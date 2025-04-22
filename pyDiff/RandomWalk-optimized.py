import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numba

'''

USARE QUESTO FILE PER PROVARE LE OTTIMIZZAZIONI:
Profile first if unsure about bottlenecks.

Implement multiprocessing for running trials in parallel (best for overall time to get averaged results).

Consider pre-computing probabilities if the disorder function is complex and called very often.

Apply Numba (@njit) to random_walk_ordered and potentially _get_grid_index_fast (if using direct calculation).
 
Try Numba on random_walk_disordered if you can refactor/precompute probabilities.

Optimize grid indexing (_get_grid_index_fast) if the grid is uniform.



Vectorize random_walk_ordered (compare performance with Numba version).


Created on April 8, 2025 by Luca Sfriso
THIS CODE SIMULATES A DISCRETE RANDOM WALK BOTH ORDER AND DISORDERED FOR A GIVEN NUMBER OF PARTICLES AND A GIVEN NUMBER OF STEPS. 
THE QUENCHED DISORDER IS SIMULATED BY CHANGING THE PROBABILITIES OF MOVING ALONG TEH COORDINATE DIRECTIONS AND RESTING IN THE SAME PLACE.


DA CAPIRE: COME USARE sys, CALCOLO DEL COEFFICIENTE DI DIFFUSIONE

UPDATE 8/04: aggiunti random_walk_ordered, ordered_trajectories, compute_msd, animate_trajectories

UPDATE 10/04: aggiornato compute_msd in modo che calcoli lungo tutte le direzioni x e y. aggiunto il disordine:
con gaussiane e rette ok, esponenziali dà problemi lungo x

UPDATE 11/04: aggiunti gli istogrammi, trovato errore nella definizione della meshgrid

UPDATE 14/04: modificata compute_D in modo che faccia un interpolazione su log log e i plot log log. NON TORNA ANCORA IL COEFFIENTE DI DIFFUSIONE CON QUELLO TEORICO

UPDATE 15/04: NON SI RIESCE A INSERIRE IL RESTING TIME IN MODO SENSATO: PLOT MSD SENZA SENSO (risolto!!)

UPDATE 16/04: Inserito il resting time. 

UPDATE 17/04: Inserite le condizioni al contorno periodiche e la funzione _check_out_of_bounds che controlla quali walker escono dalla meshgrid
                NON TORNA ANCORA IL COEFFIENTE DI DIFFUSIONE CON QUELLO TEORICO

UPDATE 18/04: Inserito il check di robustezza per capire quanto sono distante dai punti di griglia
            Inserito un modo per fare i grafici delle traiettorie sulle griglie. NON TORNA ANCORA IL COEFFIENTE DI DIFFUSIONE CON QUELLO TEORICO


UPDATE 19/04: Eseguiti test per capire se il resting time funziona. 

UPDATE 20/04: Eseguiti test per capire se il resting time funziona: aggiunta la possibilità di mediare su più run
            distinte e con resting time diversi: gaussiano, esponenziale, uniforme, plateau, multi center


TODO:
    OTTIMIZZARE E FARE NUOVI TEST
    USARE SLOPE TEST
    CAPIRE CON GIT
    CAPIRE COME USARE IL PROGRAMMA PER FARE PREDIZIONI FISICHE SUL COEFFICIENTE DI DIFFUSIONE ETC...
    NUMBA(?)

'''
@numba.njit(cache=True)
def _calculate_grid_index_fast_numba(x, y, x_min_grid, y_min_grid, dx, dy, nx, ny):
    """
    Calculates the nearest grid indices (y_idx, x_idx) for a point (x, y)
    on a uniform grid using direct calculation. Numba-optimized. Uses min/max instead of np.clip.

    Args:
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        x_min_grid (float): Minimum x-coordinate of the grid centers.
        y_min_grid (float): Minimum y-coordinate of the grid centers.
        dx (float): Grid spacing in x-direction (must be > 0).
        dy (float): Grid spacing in y-direction (must be > 0).
        nx (int): Number of grid points in x-direction.
        ny (int): Number of grid points in y-direction.

    Returns:
        tuple[int, int]: The (y_index, x_index) of the nearest grid point.
    """
    fx = (x - x_min_grid) / dx
    fy = (y - y_min_grid) / dy

    # Round to nearest integer index
    # Cast to int64 first to ensure integer type before potential clipping
    x_index = np.int64(np.round(fx))
    y_index = np.int64(np.round(fy))

    # --- Use min/max instead of np.clip ---
    # Ensure index is not less than 0
    x_index = max(0, x_index)
    y_index = max(0, y_index)

    # Ensure index is not more than nx-1 or ny-1
    x_index = min(x_index, nx - 1)
    y_index = min(y_index, ny - 1)
    # --------------------------------------

    return y_index, x_index



class RandomWalk:
    def __init__(self, interpolation=False, disorder_function=None, num_steps=1000, step=0.001, num_walkers=100,
                 xv=np.meshgrid(np.linspace(-1, 1, 2000), np.linspace(-1, 1, 2000))[0],
                 yv=np.meshgrid(np.linspace(-1, 1, 2000), np.linspace(-1, 1, 2000))[1], dt=0.0001):
        self.num_walkers = num_walkers
        self.num_steps = num_steps
        self.xv = xv
        self.yv = yv
        self.ny, self.nx = self.xv.shape
        # --- Setup for _get_grid_index_fast ---
        # Check if grid is uniform and calculate parameters if possible
        self.is_uniform_grid = False  # Flag
        if self.nx > 1 and self.ny > 1:
            self.x_coords = self.xv[0, :]
            self.y_coords = self.yv[:, 0]
            # Check for uniform spacing (within tolerance)
            dxs = np.diff(self.x_coords)
            dys = np.diff(self.y_coords)
            if np.allclose(dxs, dxs[0]) and np.allclose(dys, dys[0]):
                self.is_uniform_grid = True
                self.x_min_grid = self.x_coords[0]
                self.y_min_grid = self.y_coords[0]
                self.dx = dxs[0]
                self.dy = dys[0]
                # Ensure dx/dy are not zero
                self.dx = self.dx if abs(self.dx) > 1e-15 else 1.0
                self.dy = self.dy if abs(self.dy) > 1e-15 else 1.0
                print("Uniform grid detected. Using fast indexing.")
            else:
                print("Non-uniform grid detected. Will use robust indexing.")
        elif self.nx == 1 or self.ny == 1:
            # Handle 1D grid cases if necessary, or default to robust
            print("Grid is 1D or single point. Using robust indexing.")
            pass  # Fallback to robust indexing below
        else:
            print("Grid has zero dimensions? Using robust indexing.")
        self.box_size = np.array([xv.shape[0], yv.shape[0]])
        self.dt = dt
        self.positions = np.zeros((self.num_walkers, 2))
        self.initial_positions = np.copy(self.positions)
        # self.initial_positions=np.zeros((self.num_walkers,2))+0.5
        self.step = step
        self.all_positions = [np.copy(self.positions)]
        self.time = np.arange(self.num_steps + 1)
        self.disorder_function = disorder_function if disorder_function is not None else self._default_disorder_function
        self.x_min = np.min(self.xv)
        self.x_max = np.max(self.xv)
        self.y_min = np.min(self.yv)
        self.y_max = np.max(self.yv)
        self.out_of_bounds_walkers = set()

    def apply_pbc(self):
        """Apply periodic boundary conditions to keep particles inside the simulation box."""
        self.positions = (self.positions + self.box_size / 2) % self.box_size - self.box_size / 2

    def _get_grid_index(self, x, y):
        """
        Selects the appropriate grid index function based on uniformity.
        """
        if self.is_uniform_grid:
            return self._get_grid_index_fast(x, y)
        else:
            # Fallback to the original robust method if grid is not uniform
            # Ensure _get_grid_index_robust still exists
            return self._get_grid_index_robust(x, y)

    def _get_grid_index_robust(self, x, y):
        """Find the indices of the nearest grid point in a robust way."""

        # Calculate distances to all grid points
        x_distances = np.abs(self.xv[0, :] - x)
        y_distances = np.abs(self.yv[:, 0] - y)

        # Find the indices of the minimum distances
        x_index = np.argmin(x_distances)
        y_index = np.argmin(y_distances)

        # Robustness check (optional but recommended)
        nearest_x = self.xv[0, x_index]
        nearest_y = self.yv[y_index, 0]
        distance = np.sqrt((nearest_x - x) ** 2 + (nearest_y - y) ** 2)
        if distance > self.step * 1.1:  # Allow for small floating-point errors
            print(f"WARNING: Large distance to nearest grid point: {distance}")
            print(f"  x: {x}, y: {y}, nearest_x: {nearest_x}, nearest_y: {nearest_y}")

        return y_index, x_index

    def _get_grid_index_fast(self, x, y):
        """
        Wrapper method to call the Numba-optimized fast grid index calculation.
        Assumes a uniform grid and that grid parameters are set in __init__.
        """
        # No need for the check here if we trust the flag set in __init__
        # and only call this method when self.is_uniform_grid is True

        # Call the external Numba function, passing instance attributes
        return _calculate_grid_index_fast_numba(
            x, y,
            self.x_min_grid, self.y_min_grid,
            self.dx, self.dy,
            self.nx, self.ny
        )
    """
    #DEBUGGING
    def random_walk_disordered(self):
        print("--- random_walk_disordered ---")
        steps = np.zeros((self.num_walkers, 2))
        print(f"  Initialized steps: {steps.shape}, first value: {steps[0, 0]}")
        #moved = False

        for i in range(self.num_walkers):
            print(f"  Using disorder_function: {self.disorder_function.__name__}")
            #print(f"  --- Walker {i} ---")
            current_pos = self.positions[i]
            #print(f"    current_pos: {current_pos}")
            y_idx, x_idx = self._get_grid_index(current_pos[0], current_pos[1])
            #print(f"    y_idx: {y_idx}, x_idx: {x_idx}")
            grid_value_x = self.xv[y_idx, x_idx]
            grid_value_y = self.yv[y_idx, x_idx]
            # print(f"    grid_value_x: {grid_value_x}, grid_value_y: {grid_value_y}")
            probabilities = self.disorder_function(grid_value_x, grid_value_y)
            print(f"    probabilities: {probabilities}")

            # Ensure probabilities are valid (basic check)
            if len(probabilities) != 5 or not np.isclose(np.sum(probabilities), 1.0) or np.any(probabilities < 0):
                probabilities = self._default_disorder_function(grid_value_x, grid_value_y)
                #print(f"Warning: Invalid probabilities at ({grid_value_x:.3f}, {grid_value_y:.3f}), using default.")

            rest_prob = probabilities[4]
            print(f"    rest_prob: {rest_prob}")# Get resting probability from function



            rand_val=np.random.rand()
            #print(f"    rand_val: {rand_val}")
            if rand_val >= rest_prob:
                print("    Walker is moving")
                #moved = True
                direction_probs = probabilities[:4]
                direction_probs_sum = np.sum(direction_probs)

                if direction_probs_sum > 0:  # Avoid division by zero
                    direction_probs_normalized = direction_probs / direction_probs_sum
                else:
                    direction_probs_normalized = np.array([0.25, 0.25, 0.25, 0.25])  # Default uniform

                direction_choice = np.random.choice(4, p=direction_probs_normalized)
                # direction_choice = np.random.choice(4, p=probabilities[:4])
                #print(f"    direction_choice: {direction_choice}")

                if direction_choice == 0:  # +x
                    steps[i, 0] = self.step
                    #print(f"      Moving +x: {self.step}")
                elif direction_choice == 1:  # -x
                    steps[i, 0] = -self.step
                    # print(f"      Moving -x: {self.step}")
                elif direction_choice == 2:  # +y
                    steps[i, 1] = self.step
                    # print(f"      Moving +y: {self.step}")
                elif direction_choice == 3:  # -y
                    steps[i, 1] = -self.step
                    # print(f"      Moving -y: {self.step}")

                self.positions[i] += steps[i]

            #else: print("    Walker is resting")

            print(f"    self.positions[{i}]: {self.positions[i]}")

        print(f"  Final self.positions: {self.positions}")
        return self.positions


    """

    def random_walk_disordered(self):
        #print("--- random_walk_disordered ---")
        #print(f"  Using disorder_function: {self.disorder_function.__name__}")

        steps = np.zeros((self.num_walkers, 2))

        for i in range(self.num_walkers):

            '''
            print(f"  --- Walker {i} ---")
            current_pos = self.positions[i]
            print(f"    current_pos: {current_pos}")
            '''

            # y_idx, x_idx = self._get_grid_index(self.positions[i, 0], self.positions[i, 1])
            y_idx, x_idx = self._get_grid_index(self.positions[i, 0], self.positions[i, 1])

            # DEBUGGING
            '''
            print(f"    y_idx: {y_idx}, x_idx: {x_idx}")
            grid_value_x = self.xv[y_idx, x_idx]
            grid_value_y = self.yv[y_idx, x_idx]
            print(f"    grid_value_x: {grid_value_x}, grid_value_y: {grid_value_y}")
            '''

            probabilities = self.disorder_function(self.xv[y_idx, x_idx], self.yv[y_idx, x_idx])

            # Ensure probabilities are valid
            if len(probabilities) != 5 or not np.isclose(probabilities.sum(), 1.0) or np.any(probabilities < 0):
                probabilities = self._default_disorder_function(self.xv[y_idx, x_idx], self.yv[y_idx, x_idx])

            rest_prob = probabilities[4]

            if np.random.rand() >= rest_prob:
                direction_probs = probabilities[:4]
                direction_probs_sum = np.sum(direction_probs)

                if direction_probs_sum > 0:  # Avoid division by zero
                    direction_probs_normalized = direction_probs / direction_probs_sum
                else:
                    direction_probs_normalized = np.array([0.25, 0.25, 0.25, 0.25])  # Default uniform

                direction_choice = np.random.choice(4, p=direction_probs_normalized)
                # direction_choice = np.random.choice(4, p=probabilities[:4])

                if direction_choice == 0:  # +x
                    steps[i, 0] = self.step
                elif direction_choice == 1:  # -x
                    steps[i, 0] = -self.step
                elif direction_choice == 2:  # +y
                    steps[i, 1] = self.step
                elif direction_choice == 3:  # -y
                    steps[i, 1] = -self.step

                self.positions[i] += steps[i]  # Update position for THIS walker ONLY

            # self.apply_pbc()

        return self.positions

    def random_walk_disordered_interpolated(self):
        """
        Simulates a disordered random walk with bilinear interpolation of probabilities.
        """
        steps = np.zeros((self.num_walkers, 2))

        for i in range(self.num_walkers):
            x = self.positions[i, 0]
            y = self.positions[i, 1]

            # 1. Find the surrounding grid indices
            x_indices = np.where((self.xv[0, :] <= x))[0]
            y_indices = np.where((self.yv[:, 0] <= y))[0]

            if len(x_indices) == 0:
                x_index_left = 0
                x_index_right = 0
            elif x_indices[-1] == self.xv.shape[1] - 1:
                x_index_left = self.xv.shape[1] - 1
                x_index_right = self.xv.shape[1] - 1
            else:
                x_index_left = x_indices[-1]
                x_index_right = x_index_left + 1

            if len(y_indices) == 0:
                y_index_bottom = 0
                y_index_top = 0
            elif y_indices[-1] == self.yv.shape[0] - 1:
                y_index_bottom = self.yv.shape[0] - 1
                y_index_top = y_index_bottom + 1
            else:
                y_index_bottom = y_indices[-1]
                y_index_top = y_index_bottom + 1

            # 2. Get the coordinates of the surrounding grid points
            x1 = self.xv[0, x_index_left]
            x2 = self.xv[0, x_index_right]
            y1 = self.yv[y_index_bottom, 0]
            y2 = self.yv[y_index_top, 0]

            # 3. Get the probabilities at the surrounding grid points
            prob_Q11 = self.disorder_function(x1, y1)  # Bottom-left
            prob_Q12 = self.disorder_function(x1, y2)  # Top-left
            prob_Q21 = self.disorder_function(x2, y1)  # Bottom-right
            prob_Q22 = self.disorder_function(x2, y2)  # Top-right

            # 4. Perform bilinear interpolation
            if x1 == x2 or y1 == y2:
                interpolated_probs = prob_Q11  # Or any of the other prob_Q's
            else:
                interpolated_probs = (
                        prob_Q11 * ((x2 - x) * (y2 - y)) / ((x2 - x1) * (y2 - y1)) +
                        prob_Q21 * ((x - x1) * (y2 - y)) / ((x2 - x1) * (y2 - y1)) +
                        prob_Q12 * ((x2 - x) * (y - y1)) / ((x2 - x1) * (y2 - y1)) +
                        prob_Q22 * ((x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1))
                )

            # 5. Choose direction based on interpolated probabilities
            rest_prob = interpolated_probs[4]
            if np.random.rand() >= rest_prob:
                direction_choice = np.random.choice(4, p=interpolated_probs[:4])

                if direction_choice == 0:  # +x
                    steps[i, 0] = self.step
                elif direction_choice == 1:  # -x
                    steps[i, 0] = -self.step
                elif direction_choice == 2:  # +y
                    steps[i, 1] = self.step
                elif direction_choice == 3:  # -y
                    steps[i, 1] = -self.step

                self.positions[i] += steps[i]
        # self.apply_pbc()

        return self.positions

    def _default_disorder_function(self, x, y):
        """Default: Uniform probability distribution (no spatial disorder)."""
        return np.array([0.25, 0.25, 0.25, 0.25, 0])  # [+x, -x, +y, -y, rest]




    def random_walk_ordered(self):
        """
        Calculates the positions of the walkers at a single time step

        """
        steps = np.zeros((self.num_walkers, 2))
        directions = np.random.randint(0, 4, self.num_walkers)  # 0: +x, 1: -x, 2: +y, 3: -y

        for i in range(self.num_walkers):
            if directions[i] == 0:
                steps[i, 0] = self.step
            elif directions[i] == 1:
                steps[i, 0] = -self.step
            elif directions[i] == 2:
                steps[i, 1] = self.step
            elif directions[i] == 3:
                steps[i, 1] = -self.step

        self.positions += steps

        '''  
        ALSO A DIAGONAL MOVEMENT
        bx[b[:, 0] > 0.5] = self.step  # Move +step if random number > 0.5
        bx[b[:, 0] <= 0.5] = -self.step  # Move -step if random number <= 0.5
        by[b[:, 1] > 0.5] = self.step
        by[b[:, 1] <= 0.5] = -self.step


        self.positions[:, 0] += bx.flatten()  # Update x-coordinates for all walkers
        self.positions[:, 1] += by.flatten()  # Update y-coordinates for all walkers
        '''

        return self.positions

    def _check_out_of_bounds(self, positions):
        """
        Checks if any walkers are out of the grid and prints a warning.

        Args:
            positions (numpy.ndarray): A 2D array of walker positions (shape: (num_walkers, 2)).
        """
        out_of_bounds = np.where(
            (positions[:, 0] < self.x_min) | (positions[:, 0] > self.x_max) |
            (positions[:, 1] < self.y_min) | (positions[:, 1] > self.y_max)
        )[0]
        if len(out_of_bounds) > 0:
            print("Warning: Some walkers are out of bounds!")
            print("Out-of-bounds walker indices:", out_of_bounds)
            self.out_of_bounds_walkers.update(out_of_bounds)  # Update the set

    def trajectories(self, use_disorder=False, intepolation=False):
        self.all_positions = [np.copy(self.initial_positions)]  # Reset positions
        for _ in range(self.num_steps):
            if use_disorder:
                if intepolation:
                    self.random_walk_ordered()
                else:
                    self.random_walk_disordered()
            else:
                self.random_walk_ordered()
            # self._check_out_of_bounds(self.positions)
            self.all_positions.append(np.copy(self.positions))
        self.all_positions = np.array(self.all_positions)
        return self.all_positions

    '''
    def ordered_trajectories(self):
        """
        Calculates the positions of all the random walkers following the ordered random walk.

        """
        for _ in range(self.num_steps):
            self.random_walk_ordered()
            self.all_positions.append(np.copy(self.positions))

        self.all_positions = np.array(self.all_positions)

        return self.all_positions
    '''

    def animate_trajectory(self, walker_index=0, interval=100, save_animation=False, filename="random_walk.gif"):
        """
        Animates the trajectory of a single random walker.

        Args:
            walker_index (int): The index of the walker to animate (default: 0).
            interval (int): The delay between frames in milliseconds (default: 100).
        """
        if not hasattr(self, 'all_positions'):
            print("Error: Run the simulation first using run_simulation()")
            return

        if walker_index >= self.num_walkers:
            print(f"Error: Walker index {walker_index} is out of bounds (0 to {self.num_walkers - 1}).")
            return

        fig, ax = plt.subplots()
        ax.set_xlim(np.min(self.all_positions[:, walker_index, 0]), np.max(self.all_positions[:, walker_index, 0]))
        ax.set_ylim(np.min(self.all_positions[:, walker_index, 1]), np.max(self.all_positions[:, walker_index, 1]))
        ax.set_xlabel("X Position")
        ax.set_ylabel("Y Position")
        ax.set_title(f"Trajectory of Walker {walker_index}")
        ax.grid(visible=True, which='major', color='black', linestyle='-')
        ax.grid(visible=True, which='minor', color='black', linestyle='--')
        plt.minorticks_on()
        line, = ax.plot([], [], lw=2)
        point, = ax.plot([], [], 'ro', markersize=8)  # Current position

        def update(frame):
            xdata = self.all_positions[:frame, walker_index, 0]
            ydata = self.all_positions[:frame, walker_index, 1]
            line.set_data(xdata, ydata)
            point.set_data(xdata[-1:], ydata[-1:])  # Update current position marker
            return line, point

        ani = animation.FuncAnimation(fig, update, frames=len(self.all_positions), interval=interval, blit=True)

        if save_animation:
            try:
                ani.save(filename, writer='ffmpeg')
            except:
                print("FFmpeg not found. Saving as GIF instead.")
                ani.save("random_walk.gif", writer='pillow')
        else:
            plt.show()

    def compute_msd(self, direction='all'):
        """
        Compute the mean squared displacement over time.

        Args:
            direction (str, optional): The direction along which to calculate the MSD.
                Options are 'x', 'y', or 'all' (default).

        Returns:
            numpy.ndarray: The mean squared displacement over time.

            plottare log log il msd
        """
        all_positions_array = np.array(self.all_positions)
        displacement = all_positions_array - self.initial_positions[np.newaxis, :, :]

        if direction == 'x':
            squared_displacement = displacement[:, :, 0] ** 2  # Only x-component
        elif direction == 'y':
            squared_displacement = displacement[:, :, 1] ** 2  # Only y-component
        elif direction == 'all':
            squared_displacement = np.sum(displacement ** 2, axis=2)  # Both x and y
        else:
            raise ValueError("Invalid direction. Choose 'x', 'y', or 'all'.")

        msd = np.mean(squared_displacement, axis=1)
        return msd

    def compute_D(self):
        """
        computes the Diffusion coefficient by interpolating the msd (CAPIRE COME FARE SE MSD NON è LINEARE)
            capire se usare Dl^2/(2*Dt)
            plottare log log il msd e fittare il log log

        """
        msd = self.compute_msd(direction='all')
        time = self.time

        # Ensure that time and msd are positive for log-log plot
        valid_indices = (time > 0) & (msd > 0)
        time_valid = time[valid_indices]
        msd_valid = msd[valid_indices]
        Dp = self.step ** 2 / (2 * self.dt)  # QUELLO GIUSTO?
        D = 0.25 * msd_valid / time_valid

        if np.sum(valid_indices) < 2:
            print("Warning: Not enough valid data points for log-log plot.")
            return

        # Perform the linear fit in log-log space
        log_time = np.log(time_valid)
        log_msd = np.log(msd_valid)
        try:
            p = np.polyfit(log_time, log_msd, 1)
            slope, intercept = p
        except np.RankWarning:
            print("Warning: Log-log fit may be poorly conditioned.")
            slope, intercept = 1.0, 0.0  # Default values
        except np.linalg.LinAlgError:
            print("Warning: Linear algebra error during log-log fit.")
            return

        '''
        D_All=np.polyfit(self.time,self.compute_msd(),1)[0]*0.25
        D_x=np.polyfit(self.time,self.compute_msd(direction='x'),1)[0]*0.5
        D_y = np.polyfit(self.time, self.compute_msd(direction='y'), 1)[0] * 0.5
        '''
        return Dp, D, slope, intercept

    def plot_position_histograms(self, time_step, walker_indices=None):
        """
        Plots histograms of walker positions at a given time step.

        Args:
            time_step (int): The time step for which to plot the histograms.
            walker_indices (list, optional): A list of walker indices to include in the histograms.
                                           If None, all walkers are included (default: None).
        """

        if not hasattr(self, 'all_positions'):
            raise ValueError("Run the simulation first using ordered_trajectories()")

        if time_step < 0 or time_step >= self.all_positions.shape[0]:
            raise ValueError(f"Invalid time step: {time_step}. Must be between 0 and {self.all_positions.shape[0] - 1}")

        positions_at_time = self.all_positions[time_step]  # (num_walkers, 2) array

        if walker_indices is None:
            x_positions = positions_at_time[:, 0]
            y_positions = positions_at_time[:, 1]
        else:
            x_positions = positions_at_time[walker_indices, 0]
            y_positions = positions_at_time[walker_indices, 1]

        fig, axs = plt.subplots(1, 2, figsize=(10, 5))

        axs[0].hist(x_positions, bins=10, color='skyblue', edgecolor='black')
        axs[0].set_xlabel("X Position")
        axs[0].set_ylabel("Frequency")
        axs[0].set_title(f"X Positions at Time Step {time_step}")

        axs[1].hist(y_positions, bins=10, color='lightgreen', edgecolor='black')
        axs[1].set_xlabel("Y Position")
        axs[1].set_ylabel("Frequency")
        axs[1].set_title(f"Y Positions at Time Step {time_step}")

        plt.tight_layout()
        plt.show()

    def plot_trajectory_on_grid_zoomed(self, walker_index=0, num_steps_to_plot=None, zoom_factor=1):
        """
        Plots the trajectory of a single walker, zooming in on the path.

        Args:
            walker_index (int, optional): The index of the walker to plot (default: 0).
            num_steps_to_plot (int, optional): The number of steps to plot. If None, plots the entire trajectory (default: None).
            zoom_factor (float, optional): A factor to control the zoom level (default: 1.2).
                                         Values > 1 zoom out, values < 1 zoom in further.
        """

        if not hasattr(self, 'all_positions'):
            raise ValueError("Run the simulation first using trajectories()")

        if walker_index >= self.num_walkers:
            raise ValueError(f"Walker index {walker_index} out of range (0 to {self.num_walkers - 1})")

        positions = self.all_positions[:, walker_index, :]  # Get trajectory of the specified walker

        if num_steps_to_plot is None:
            plot_positions = positions
        else:
            plot_positions = positions[:num_steps_to_plot]

        # Calculate plot limits based on walker's trajectory
        x_min = np.min(plot_positions[:, 0])
        x_max = np.max(plot_positions[:, 0])
        y_min = np.min(plot_positions[:, 1])
        y_max = np.max(plot_positions[:, 1])

        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        x_range = (x_max - x_min) * zoom_factor / 2
        y_range = (y_max - y_min) * zoom_factor / 2

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlim(x_center - x_range, x_center + x_range)
        ax.set_ylim(y_center - y_range, y_center + y_range)

        # Plot the grid (only within the zoomed view)
        x_grid_min = x_center - x_range
        x_grid_max = x_center + x_range
        y_grid_min = y_center - y_range
        y_grid_max = y_center + y_range

        x_grid_indices = np.where((self.xv[0, :] >= x_grid_min) & (self.xv[0, :] <= x_grid_max))[0]
        y_grid_indices = np.where((self.yv[:, 0] >= y_grid_min) & (self.yv[:, 0] <= y_grid_max))[0]

        ax.plot(self.xv[y_grid_indices, :][:, x_grid_indices], self.yv[y_grid_indices, :][:, x_grid_indices], 'k-',
                linewidth=0.5, alpha=0.5)  # Vertical
        ax.plot(self.xv[y_grid_indices, :][:, x_grid_indices].T, self.yv[y_grid_indices, :][:, x_grid_indices].T, 'k-',
                linewidth=0.5, alpha=0.5)  # Horizontal

        # Plot the walker's trajectory
        ax.plot(plot_positions[:, 0], plot_positions[:, 1], 'r-', label=f"Walker {walker_index} Trajectory")
        ax.plot(plot_positions[0, 0], plot_positions[0, 1], 'go', markersize=8, label="Start")
        ax.plot(plot_positions[-1, 0], plot_positions[-1, 1], 'mo', markersize=8, label="End")

        ax.set_xlabel("X Position")
        ax.set_ylabel("Y Position")
        ax.set_title(f"Random Walk Trajectory (Zoomed) (Walker {walker_index})")
        ax.legend()
        ax.set_aspect('equal', 'box')
        plt.show()


def my_spatial_disorder(x, y):
    """LAPLACE
    if(x<0):
        prob_plus_x=0.05*np.exp(0.1*x)
        prob_minus_x =0.05*np.exp(0.1*x)
    else:
        prob_plus_x=0.05*np.exp(-0.1*x)
        prob_minus_x=0.05*np.exp(-0.1*x)

    if (y < 0):
        prob_plus_y = 0.25 * np.exp(0.5 * y)
        prob_minus_y = 0.25 * np.exp(0.5 * y)
    else:
        prob_plus_y= 0.25 * np.exp(-0.5 * y)
        prob_minus_y = 0.25* np.exp(-0.5 * y)
    """
    '''
    #GAUSSIANA
    prob_plus_x=1/np.sqrt(np.pi*(2*0.1))*np.exp(-x**2/(2*0.1))
    prob_minus_x = 1 / np.sqrt(np.pi * (2*0.1)) * np.exp(-x ** 2 / (2*0.1))
    prob_plus_y = 1 / np.sqrt(np.pi * (2*0.1)) * np.exp(-y ** 2 / (2*0.1))
    prob_minus_y = 1 / np.sqrt(np.pi * (2*0.1)) * np.exp(-y ** 2 / (2*0.1))
    rest_prob =1 / np.sqrt(np.pi * (4 * 0.5)) * np.exp(-(y ** 2+x**2) / (4 * 0.5))

    rest_prob =(1/(2*np.pi*0.1*0.5)*np.exp(-(x**2)/(2*0.1**2)-y**2/(2*0.5**2)))*0.2 #0.2 serve per normalizzarle tutte
    prob_plus_x=1/np.sqrt(np.pi*(2*0.1))*np.exp(-x**2/(2*0.1**2))*0.2
    prob_minus_x = 1 / np.sqrt(np.pi * (2*0.1)) * np.exp(-x ** 2 / (2*0.1**2))*0.2
    prob_plus_y = 1 / np.sqrt(np.pi * (2*0.5)) * np.exp(-y ** 2 / (2*0.5**2))*0.2
    prob_minus_y = 1 / np.sqrt(np.pi * (2*0.5)) * np.exp(-y ** 2 / (2*0.5**2))*0.2
    '''
    #
    '''

    '''
    rest_prob = 0.9
    prob_plus_x = 0.25 * (1 - rest_prob)
    prob_minus_x = 0.25 * (1 - rest_prob)
    prob_plus_y = 0.25 * (1 - rest_prob)
    prob_minus_y = 0.25 * (1 - rest_prob)

    # probs = np.clip([prob_plus_x, prob_minus_x, prob_plus_y, prob_minus_y,rest_prob], 0, 1)
    probs = [prob_plus_x, prob_minus_x, prob_plus_y, prob_minus_y, rest_prob]

    if np.sum(probs) > 0:
        probs = probs / np.sum(probs)
    else:
        probs = [0.25, 0.25, 0.25, 0.25, 0]

    return probs


def gaussian_rest_prob(x, y, sigma=0.1):
    """Gaussian resting probability centered at the origin."""
    return np.array([0.2, 0.2, 0.2, 0.2, np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2))])


def gaussian_rest_prob_streght(x, y, sigma=0.1, rest_strength=0.95):
    """Gaussian resting probability centered at the origin."""
    rest_prob = rest_strength * np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2))
    return np.clip([0.2, 0.2, 0.2, 0.2, rest_prob], 0, 1)


def uniform_rest_prob(x, y, rest_level=0.9):
    """Uniform resting probability across the grid."""
    return np.array([0.25 * (1 - rest_level), 0.25 * (1 - rest_level),
                     0.25 * (1 - rest_level), 0.25 * (1 - rest_level),
                     rest_level])


def plateau_rest_prob(x, y, plateau_radius=0.001, max_rest=0.5, decay_rate=5.0):
    """Resting probability has a plateau near the origin."""
    distance_from_origin = np.sqrt(x ** 2 + y ** 2)
    if distance_from_origin <= plateau_radius:
        rest_prob = max_rest
    else:
        rest_prob = max_rest * np.exp(-decay_rate * (distance_from_origin - plateau_radius))
    base_prob = (1 - rest_prob) / 4
    return np.clip([base_prob, base_prob, base_prob, base_prob, rest_prob], 0, 1)


def multi_center_rest_prob(x, y):
    """High resting probability around multiple centers."""
    center1_rest = 0.5 * np.exp(-5 * ((x - 0.001) ** 2 + (y - 0.001) ** 2))
    center2_rest = 0.5 * np.exp(-5 * ((x + 0.001) ** 2 + (y + 0.001) ** 2))
    rest_prob = np.clip(center1_rest + center2_rest, 0, 0.9)
    base_prob = (1 - rest_prob) / 4
    return np.clip([base_prob, base_prob, base_prob, base_prob, rest_prob], 0, 1)


def exponential_rest_prob(x, y, decay_rate=0.1, max_rest=0.95):
    """Resting probability decreases exponentially from the origin."""
    rest_prob = max_rest * np.exp(-decay_rate * np.sqrt(x ** 2 + y ** 2))
    return np.clip([0.2 * (1 - max_rest), 0.2 * (1 - max_rest),
                    0.3 * (1 - max_rest), 0.3 * (1 - max_rest),
                    rest_prob], 0, 1)


def boundary_dependent_rest_prob(x, y, boundary_strength=0.5):
    """Resting probability increases near the boundaries. QUESTO NON HA TANTO SENSO"""
    rest_prob_x = boundary_strength * (np.exp(-5 * (1 - np.abs(x))) + np.exp(-5 * (1 + x)))
    rest_prob_y = boundary_strength * (np.exp(-5 * (1 - np.abs(y))) + np.exp(-5 * (1 + y)))
    rest_prob = np.clip(rest_prob_x + rest_prob_y, 0, 0.9)  # Combine and clip
    base_prob = (1 - rest_prob) / 4
    return np.array([base_prob, base_prob, base_prob, base_prob, rest_prob])

def corrected_gaussian_rest_prob(x, y, sigma=1, max_rest_strength=0.95):
    """
    Corrected Gaussian resting probability centered at the origin.
    Ensures the returned probabilities always sum to 1.0.

    Args:
        x (float): x-coordinate.
        y (float): y-coordinate.
        sigma (float): Standard deviation of the Gaussian distribution.
        max_rest_strength (float): The maximum resting probability at the origin (must be <= 1).

    Returns:
        numpy.ndarray: Array of 5 probabilities [p(+x), p(-x), p(+y), p(-y), p(rest)]
                       that sums to 1.0.
    """
    # Ensure max_rest_strength is valid
    max_rest_strength = min(max_rest_strength, 1.0) # Cannot be more than 1

    # Calculate the resting probability based on Gaussian decay
    rest_prob = max_rest_strength * np.exp(-(x**2 + y**2) / (2 * sigma**2))

    # Ensure rest_prob is strictly within [0, 1] after calculation
    # (Gaussian is always non-negative, clip ensures it doesn't exceed 1 if max_rest_strength > 1 was passed somehow)
    rest_prob = np.clip(rest_prob, 0.0, 1.0)

    # Calculate the total probability available for movement
    move_prob_total = 1.0 - rest_prob

    # Distribute movement probability equally among the 4 directions (+x, -x, +y, -y)
    # Handle potential floating point inaccuracies where move_prob_total might be slightly < 0
    move_prob_each = max(0.0, move_prob_total / 4.0)

    # Create the final probability array
    probs = np.array([move_prob_each, move_prob_each, move_prob_each, move_prob_each, rest_prob])

    # --- Optional: Check for debugging ---
    # if not np.isclose(np.sum(probs), 1.0):
    #     print(f"Warning: Probabilities do not sum to 1 at ({x:.3f}, {y:.3f}). Probs: {probs}, Sum: {np.sum(probs)}")
    #     # Attempt to re-normalize as a fallback, though it shouldn't be needed with this logic
    #     if np.sum(probs) > 1e-9: # Avoid division by zero
    #        probs = probs / np.sum(probs)
    #     else: # If sum is zero, default to equal probability (though rest_prob would be 1 here)
    #        probs = np.array([0.0, 0.0, 0.0, 0.0, 1.0])
    # -------------------------------------

    return probs













def main():
    num_steps = 1000
    num_trials = 1  # Set the number of independent trials
    time = np.arange(num_steps + 1)

    # Store MSD results for each trial
    msd_no_rest_trials = []
    #msd_gaussian_rest_trials = []
    # msd_high_rest_trials = []
    # msd_uniform_rest_trials = []
    # msd_exponential_rest_trials = []
    #msd_plateau_rest_trials = []
    #msd_multi_center_rest_trials = []
    msd_gaussian_rest_trials_new = []

    for _ in range(num_trials):
        # Simulate with different resting probabilities

        rw_no_rest = RandomWalk(disorder_function=None, num_steps=num_steps)
        rw_no_rest.trajectories(use_disorder=False)
        msd_no_rest_trials.append(rw_no_rest.compute_msd())
        '''
        rw_gaussian_rest = RandomWalk(disorder_function=gaussian_rest_prob_streght, num_steps=num_steps)
        rw_gaussian_rest.trajectories(use_disorder=True)
        msd_gaussian_rest_trials.append(rw_gaussian_rest.compute_msd())
        '''
        rw_gaussian_rest_new = RandomWalk(disorder_function=corrected_gaussian_rest_prob, num_steps=num_steps)
        rw_gaussian_rest_new.trajectories(use_disorder=True)
        msd_gaussian_rest_trials_new.append(rw_gaussian_rest_new.compute_msd())
        '''
        rw_high_rest = RandomWalk(disorder_function=my_spatial_disorder, num_steps=num_steps)
        rw_high_rest.trajectories(use_disorder=True)
        msd_high_rest_trials.append(rw_high_rest.compute_msd())  # Exclude time 0

        rw_uniform_rest = RandomWalk(disorder_function=uniform_rest_prob, num_steps=num_steps)
        rw_uniform_rest.trajectories(use_disorder=True)
        msd_uniform_rest_trials.append(rw_uniform_rest.compute_msd())

        rw_exponential_rest = RandomWalk(disorder_function=exponential_rest_prob, num_steps=num_steps)
        rw_exponential_rest.trajectories(use_disorder=True)
        msd_exponential_rest_trials.append(rw_exponential_rest.compute_msd())
        

        rw_plateu = RandomWalk(disorder_function=plateau_rest_prob, num_steps=num_steps)
        rw_plateu.trajectories(use_disorder=True)
        msd_plateau_rest_trials.append(rw_plateu.compute_msd())

        rw_multi_center = RandomWalk(disorder_function=multi_center_rest_prob, num_steps=num_steps)
        rw_multi_center.trajectories(use_disorder=True)
        msd_multi_center_rest_trials.append(rw_multi_center.compute_msd())
        '''
        # Calculate the average MSD over all trials
    avg_msd_no_rest = np.mean(msd_no_rest_trials, axis=0)
    #avg_msd_gaussian_rest = np.mean(msd_gaussian_rest_trials, axis=0)
    avg_msd_gaussian_rest_new = np.mean(msd_gaussian_rest_trials_new, axis=0)
    # avg_msd_high_rest = np.mean(msd_high_rest_trials, axis=0)
    # avg_msd_uniform_rest = np.mean(msd_uniform_rest_trials, axis=0)
    # avg_msd_exponential_rest = np.mean(msd_exponential_rest_trials, axis=0)
    #avg_msd_plateau = np.mean(msd_plateau_rest_trials, axis=0)
    #avg_msd_multi_center = np.mean(msd_multi_center_rest_trials, axis=0)

    '''
    # Plotting averaged MSD on a log-log scale
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(time, avg_msd_no_rest, label='No Resting (rest_prob=0)')
   # ax.loglog(time, avg_msd_gaussian_rest, label='Gaussian Resting Probability')
    ax.loglog(time, avg_msd_gaussian_rest_new, label='Gaussian Resting Probability Corrected')
    # ax.loglog(time, avg_msd_high_rest, label='High Resting Probability (rest_prob=0.9)')
    # ax.loglog(time, avg_msd_uniform_rest, label='Uniform Resting Probability')
    # ax.loglog(time, avg_msd_exponential_rest, label='Exponential Resting Probability')
   # ax.loglog(time, avg_msd_plateau, label='Plateau Resting Probability')
    #ax.loglog(time, avg_msd_multi_center, label='Multi Center Resting Probability')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('MSD')
    ax.set_title(f'Mean Squared Displacement (Log-Log Scale) - Averaged over {num_trials} Trials')
    ax.grid(True)
    ax.legend()
    plt.show()

    # Plotting averaged MSD/Time
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time[1:], avg_msd_no_rest[1:] / time[1:], label='No Resting (rest_prob=0)')
    #ax.plot(time[1:], avg_msd_gaussian_rest[1:] / time[1:], label='Gaussian Resting Probability')
    ax.plot(time[1:], avg_msd_gaussian_rest_new[1:] / time[1:], label='Gaussian Resting Probability Corrected')
    # ax.plot(time[1:], avg_msd_high_rest[1:] / time[1:], label='High Resting Probability (rest_prob=0.9)')
    # ax.plot(time[1:], avg_msd_uniform_rest[1:] / time[1:], label='Uniform Resting Probability')
    # ax.plot(time[1:], avg_msd_exponential_rest[1:] / time[1:], label='Exponential Resting Probability')
    #ax.plot(time[1:], avg_msd_plateau[1:] / time[1:], label='Plateau Resting Probability')
    #ax.plot(time[1:], avg_msd_multi_center[1:] / time[1:], label='Multi Center Resting Probability')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('MSD / Time')
    ax.set_title(f'Effective Diffusion Coefficient (MSD / Time) - Averaged over {num_trials} Trials')
    ax.grid(True)
    ax.legend()
    plt.show()
    '''
    '''
    #TEST PER IL RESTING TIME

        # Test with no resting probability
    num_steps=100
    time = np.arange(1, num_steps + 1)
    rw_no_rest = RandomWalk(disorder_function=None)
    rw_no_rest.trajectories(use_disorder=False)
    msd_no_rest = rw_no_rest.compute_msd()[1:]  # Exclude time 0

    # Test with Gaussian resting probability
    rw_gaussian_rest = RandomWalk(disorder_function=gaussian_rest_prob_streght, num_steps=num_steps)
    rw_gaussian_rest.trajectories(use_disorder=True)
    msd_gaussian_rest = rw_gaussian_rest.compute_msd()[1:]  # Exclude time 0

    # Test with very high resting probability
    rw_high_rest = RandomWalk(disorder_function=my_spatial_disorder, num_steps=num_steps)
    rw_high_rest.trajectories(use_disorder=True)
    msd_high_rest = rw_high_rest.compute_msd()[1:]  # Exclude time 0

    #Test with uniform rest prob
    rw_uniform = RandomWalk(disorder_function=uniform_rest_prob, num_steps=num_steps)
    rw_uniform.trajectories(use_disorder=True)
    msd_uniform = rw_uniform.compute_msd()[1:]  # Exclude time 0

    # Test with uniform rest prob
    rw_expo = RandomWalk(disorder_function=exponential_rest_prob, num_steps=num_steps)
    rw_expo.trajectories(use_disorder=True)
    msd_expo = rw_expo.compute_msd()[1:]  # Exclude time 0




    # Plotting MSD/Time
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time, msd_no_rest / time, label='No Resting (rest_prob=0)')
    ax.plot(time, msd_gaussian_rest / time, label='Gaussian Resting Probability')
    ax.plot(time, msd_high_rest / time, label='High Resting Probability (rest_prob=0.9)')
    ax.plot(time, msd_uniform / time, label='Uniform resting probability (rest_prob=0.9)')
    ax.plot(time, msd_expo / time, label='Exponential resting probability')


    ax.set_xlabel('Time Step')
    ax.set_ylabel('MSD / Time')
    ax.set_title('Effective Diffusion Coefficient (MSD / Time)')
    ax.grid(True)
    ax.legend()
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(time, msd_no_rest, label='No Resting (rest_prob=0)')
    ax.loglog(time, msd_gaussian_rest, label='Gaussian Resting Probability')
    ax.loglog(time, msd_high_rest, label='High Resting Probability (rest_prob=0.9)')
    ax.loglog(time, msd_uniform, label='Uniform Resting Probability')
    ax.loglog(time, msd_expo, label='Exponential Resting Probability')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('MSD')
    ax.set_title('Mean Squared Displacement (Log-Log Scale)')
    ax.grid(True)
    ax.legend()
    plt.show()
    '''

    '''
    #rw = RandomWalk(disorder_function=my_spatial_disorder)
    #rw=RandomWalk()
    rw=RandomWalk(disorder_function=gaussian_rest_prob)
    #rwo=RandomWalk(disorder_function=None)
    num_steps=rw.num_steps
    time = rw.time


    rw.trajectories(use_disorder=True)
    #rwo.trajectories(use_disorder=False)



    #rw.plot_position_histograms(time_step=100)

    msd = rw.compute_msd()
    #msdo=rwo.compute_msd()



    D=rw.compute_D()
    print("Dp=",D[0])
    print("slope (no disorder)[1]:",D[2])
    print("intercept (no disorder)[]:", np.exp(D[3])/4)

    #rw.plot_trajectory_on_grid_zoomed()



    # Plotting MSD with slope 1 line
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(time, msd, 'b-', label='MSD')

    # Generate a straight line with slope 1 in log-log space
    # y = C * t^1  => log(y) = log(C) + 1 * log(t)
    # We'll choose a C that roughly matches the MSD at later times for visual comparison
    # Find a reasonable starting point for the line (e.g., halfway through the simulation)
    mid_time_index = num_steps // 2
    C = msd[mid_time_index] / time[mid_time_index] if time[mid_time_index] > 0 else 1e-6
    slope_1_line = C * time

    ax.loglog(time, slope_1_line, 'g--', label='Slope = 1')

    ax.set_xlabel('Time Step')
    ax.set_ylabel('MSD')
    ax.set_title('Mean Squared Displacement with Slope 1 Line')
    ax.grid(True)
    ax.legend()
    plt.show()

    # Plotting MSD normalized by time
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(time[1:], msd[1:] / time[1:], 'b-', label='MSD/Time')
    #ax.plot(time[1:], msdo[1:] / time[1:], 'r-', label='MSDO/Time')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('MSD/Time')
    ax.set_title('MSD normalized by time')
    ax.grid(True)
    ax.legend()
    plt.show()


    # Plot position histograms at different times
    rw.plot_position_histograms(time_step=0)
    rw.plot_position_histograms(time_step=num_steps // 4)
    rw.plot_position_histograms(time_step=num_steps // 2)
    rw.plot_position_histograms(time_step=num_steps - 1)
    '''

    '''

    # Plot the trajectory of the first walker

    total_msd=rw.compute_msd()
    x_msd=rw.compute_msd(direction='x')
    y_msd=rw.compute_msd(direction='y')
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].minorticks_on()
    ax[0].plot(np.array(rw.all_positions)[:, 0, 0], np.array(rw.all_positions)[:, 0, 1], c='r')
    ax[0].set_xlabel("X Position")
    ax[0].set_ylabel("Y Position")
    ax[0].set_title("Trajectory of the First Walker")
    ax[0].grid(visible=True, which='major', color='black', linestyle='-')
    ax[0].grid(visible=True, which='minor', color='black', linestyle='--')

    ax[0].set_aspect('equal', adjustable='box')

    # Plot the Mean Squared Displacement over time
    ax[1].minorticks_on()
    ax[1].plot(time, total_msd, c='b')
    ax[1].set_title("Mean Squared Displacement")
    ax[1].grid(visible=True, which='major', color='black', linestyle='-')
    ax[1].grid(visible=True, which='minor', color='black', linestyle='--')
    ax[1].set_yscale('log')
    ax[1].set_xscale('log')
    plt.tight_layout()
    plt.plot(time, total_msd, label='Total MSD')
    plt.plot(time, x_msd, label='X MSD')
    plt.plot(time, y_msd, label='Y MSD')
    plt.xlabel('Time Step')
    plt.ylabel('MSD')
    plt.legend()
    plt.show()


    fig, ax = plt.subplots()
    valid_indices = (time > 0)
    time_valid = time[valid_indices]
    ax.plot(time_valid,D[1])
    plt.show()
    '''
    # rw.animate_trajectory(walker_index=0, interval=50)


if __name__ == "__main__":
    main()