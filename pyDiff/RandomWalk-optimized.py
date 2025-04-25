import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numba
import multiprocessing # Import the module
import time as timer # To time the execution
import os # To potentially get CPU count
import argparse # Import argparse
'''


Created on April 8, 2025 by Luca Sfriso
THIS CODE SIMULATES A DISCRETE RANDOM WALK BOTH ORDER AND DISORDERED FOR A GIVEN NUMBER OF PARTICLES AND A GIVEN NUMBER OF STEPS. 
THE QUENCHED DISORDER IS SIMULATED BY CHANGING THE PROBABILITIES OF MOVING ALONG TEH COORDINATE DIRECTIONS AND RESTING IN THE SAME PLACE.

USAGE EXAMPLES:

python RandomWalk-optimized.py --steps 5000 --trials 20
SIMULATES 5000 STEPS FOR 20 TRIALS WITHOUT DISORDER

python RandomWalk-optimized.py --steps 5000 --trials 20 --disorder gaussian --sigma 0.5 --max_rest 0.9
SIMULATES 5000 STEPS FOR 20 TRIALS WITH GAUSSIAN RESTING TIME DISORDER WITH SIGMA 0.5 AND MAX REST 0.9

python RandomWalk-optimized.py --steps 10000 --trials 10 --disorder uniform --rest_level 0.5 --ctrw --alpha 0.6
SIMULATES 10000 STEPS FOR 10 TRIALS WITH UNIFORM RESTING TIME DISORDER WITH REST LEVEL 0.5 AND USING THE CONTINUOS TIME RANDOM WALK WITH EXPONENT 0.6

python RandomWalk-optimized.py --steps 2000 --trials 5 --disorder plateau --plateau_radius 0.1 --max_rest 0.8 --decay_rate 10 --animate --anim_steps 1000 --anim_file plateau_walk.gif

Run simulation and show histograms for steps 1000, 5000, and 10000
python RandomWalk-optimized.py --steps 10000 --trials 10 --histograms --hist_steps 1000 5000 10000

Run Gaussian disorder and show histograms
python RandomWalk-optimized.py --steps 5000 --disorder gaussian --sigma 0.2 --max_rest 0.9 --histograms --hist_steps 500 4000

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

UPDATE 23/04: Il resting time funziona: inserito il waiting time (se il walker decide di fermarsi estrae un numero di timestep
                durante i quali rimane fermo da una powerlaw tau^-(1+alpha). Inseriti test quantitativi per verificare la subdiffusione)

TODO:
    CAPIRE COME USARE IL PROGRAMMA PER FARE PREDIZIONI FISICHE SUL COEFFICIENTE DI DIFFUSIONE ETC...

'''
@numba.njit(cache=True)
def draw_power_law_wait_time(alpha):
    v = 1.0 - np.random.rand()
    wait_steps = np.int64(np.ceil(v**(-1.0 / alpha)))
    return max(np.int64(1), wait_steps)

@numba.njit(cache=True, fastmath=True) # Consider Numba for this too
def random_walk_ordered_numba(positions, step_size, num_walkers):
    steps = np.zeros_like(positions)
    # Generate random integers 0, 1, 2, 3 for all walkers at once
    direction_choices = np.random.randint(0, 4, num_walkers)

    # Create masks for each direction
    mask_px = (direction_choices == 0)
    mask_mx = (direction_choices == 1)
    mask_py = (direction_choices == 2)
    mask_my = (direction_choices == 3)

    # Apply steps based on masks
    steps[mask_px, 0] = step_size
    steps[mask_mx, 0] = -step_size
    steps[mask_py, 1] = step_size
    steps[mask_my, 1] = -step_size

    return steps # Return only the steps
# --- Numba Kernel for STANDARD Disordered Step (Rest = 1 step) ---
@numba.njit(cache=True)
def _run_standard_disordered_step_numba(
    positions,           # Current positions (num_walkers, 2)
    step_size,           # float
    precomputed_probs,   # (ny, nx, 5) float32 array
    x_min_grid, y_min_grid, dx, dy, nx, ny # Grid parameters
    ):
    """ Numba kernel for one STANDARD disordered step. """
    num_walkers = positions.shape[0]
    steps = np.zeros_like(positions)

    for i in range(num_walkers):
        current_x = positions[i, 0]
        current_y = positions[i, 1]
        y_idx, x_idx = _calculate_grid_index_fast_numba(
            current_x, current_y,
            x_min_grid, y_min_grid, dx, dy, nx, ny)

        probabilities = precomputed_probs[y_idx, x_idx, :]
        rest_prob = probabilities[4]

        if np.random.rand() >= rest_prob: # Check if NOT resting
            direction_probs = probabilities[:4]
            direction_probs_sum = np.sum(direction_probs)
            if direction_probs_sum > 1e-9:
                choice_rand = np.random.rand() * direction_probs_sum
                p_plus_x = direction_probs[0]
                p_minus_x = direction_probs[1]
                p_plus_y = direction_probs[2]

                if choice_rand < p_plus_x: direction_choice = 0
                elif choice_rand < (p_plus_x + p_minus_x): direction_choice = 1
                elif choice_rand < (p_plus_x + p_minus_x + p_plus_y): direction_choice = 2
                else: direction_choice = 3

                if direction_choice == 0: steps[i, 0] = step_size
                elif direction_choice == 1: steps[i, 0] = -step_size
                elif direction_choice == 2: steps[i, 1] = step_size
                elif direction_choice == 3: steps[i, 1] = -step_size
        # else: If resting, steps[i,:] remains [0, 0]

    return steps


@numba.njit(cache=True)
def _run_ctrw_disordered_step_numba(
    positions,           # Input: Current (x,y) for all walkers
    wait_times,          # Input/Output: Steps remaining to wait for each walker
    step_size,           # Input: How far a walker moves if it moves
    precomputed_probs,   # Input: The table of probabilities for each grid site
    x_min_grid, y_min_grid, dx, dy, nx, ny, # Input: Grid info for indexing
    alpha                # Input: CTRW exponent for power-law waits
    ):
    """ Numba kernel for one CTRW step with power-law waits. """

    num_walkers = positions.shape[0] # How many walkers are there?

    # Initialize array to store the step [dx, dy] taken by each walker IN THIS TIMESTEP.
    # If a walker waits or fails to move, its step remains [0, 0].
    steps = np.zeros_like(positions)

    # IMPORTANT: The 'wait_times' array passed in will be modified directly by this function.

    # Loop through each walker individually
    for i in range(num_walkers):

        # --- 1. Check if the walker is currently waiting ---
        if wait_times[i] > 0:
            # If wait_times[i] is greater than 0, this walker was previously told
            # to wait and still has time left on its counter.
            wait_times[i] -= 1 # Decrement the remaining wait time by 1 step.
            # The walker does nothing else this step; its step remains [0, 0].

        # --- 2. If the walker is NOT waiting ---
        else:
            # wait_times[i] is 0, so the walker is free to act.
            # It will either attempt to move or decide to start a new rest period.

            # --- 2a. Find where the walker is on the grid ---
            current_x = positions[i, 0]
            current_y = positions[i, 1]
            # Get the grid indices (y_idx, x_idx) corresponding to the position.
            y_idx, x_idx = _calculate_grid_index_fast_numba(
                current_x, current_y, x_min_grid, y_min_grid, dx, dy, nx, ny)

            # --- 2b. Look up the probabilities for that grid site ---
            # Retrieve the 5 probabilities [p+x, p-x, p+y, p-y, p_rest]
            # from the precomputed table for this specific grid location.
            probabilities = precomputed_probs[y_idx, x_idx, :]
            rest_prob = probabilities[4] # Extract the resting probability

            # --- 2c. Decide: Rest or Move? ---
            # Generate a random float between 0.0 and 1.0
            if np.random.rand() < rest_prob:
                # --- Walker CHOOSES TO REST ---
                # Draw a waiting time 'tau' (integer >= 1) from the power-law distribution
                tau = draw_power_law_wait_time(alpha)
                # Set the walker's wait counter. The total wait is 'tau' steps.
                # Since this current step is the first step of waiting, the
                # REMAINING wait time to set on the counter is tau - 1.
                # Ensure it's not negative if tau happened to be 1.
                wait_times[i] = max(np.int64(0), tau - 1)
                # Walker doesn't move this step, steps[i,:] remains [0, 0].

            else:
                # --- Walker CHOOSES TO MOVE ---
                # Get the probabilities for the 4 move directions
                direction_probs = probabilities[:4]
                direction_probs_sum = np.sum(direction_probs) # Should equal 1.0 - rest_prob

                # Check if movement is actually possible (i.e., rest_prob wasn't 1.0)
                if direction_probs_sum > 1e-9: # Use tolerance for float comparison
                    # Choose a direction based on the relative probabilities p(+x), p(-x), p(+y), p(-y)
                    # This implements np.random.choice(4, p=normalized_probs) efficiently for Numba:
                    choice_rand = np.random.rand() * direction_probs_sum # Random number scaled to total move prob
                    p_plus_x  = direction_probs[0]
                    p_minus_x = direction_probs[1]
                    p_plus_y  = direction_probs[2]

                    # Check cumulatively which direction bin the random number falls into
                    if choice_rand < p_plus_x:
                        direction_choice = 0 # +x
                    elif choice_rand < (p_plus_x + p_minus_x):
                        direction_choice = 1 # -x
                    elif choice_rand < (p_plus_x + p_minus_x + p_plus_y):
                        direction_choice = 2 # +y
                    else:
                        direction_choice = 3 # -y

                    # Assign the actual step [dx, dy] based on the chosen direction
                    if direction_choice == 0: steps[i, 0] = step_size
                    elif direction_choice == 1: steps[i, 0] = -step_size
                    elif direction_choice == 2: steps[i, 1] = step_size
                    elif direction_choice == 3: steps[i, 1] = -step_size

                # If direction_probs_sum was ~0 (rest_prob was ~1), the walker cannot move.
                # steps[i,:] remains [0, 0].
                # In either move case (moved or couldn't move because rest_prob=1),
                # the wait_times[i] remains 0 because the walker didn't *choose* to rest.

    # After looping through all walkers, return the array of steps taken in this dt
    return steps


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
    def __init__(self, interpolation=False,store_history=False, disorder_function=None, num_steps=1000, step=0.001, num_walkers=100,
                 xv=np.meshgrid(np.linspace(-1, 1, 2000), np.linspace(-1, 1, 2000))[0],
                 yv=np.meshgrid(np.linspace(-1, 1, 2000), np.linspace(-1, 1, 2000))[1], dt=0.0001,disorder_params={'sigma': 1.0, 'max_rest_strength': 0.95},use_ctrw=False):
        self.num_walkers = num_walkers
        self.num_steps = num_steps
        # --- Grid Setup ---
        if xv is None or yv is None:
            # Default grid if none provided
            print("Warning: No grid provided, using default 100x100 grid.")
            xv, yv = np.meshgrid(np.linspace(-1, 1, 100), np.linspace(-1, 1, 100))
        self.xv = xv
        self.yv = yv
        self.ny, self.nx = self.xv.shape
        self.x_min = np.min(self.xv)
        self.x_max = np.max(self.xv)
        self.y_min = np.min(self.yv)
        self.y_max = np.max(self.yv)
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
        self.disorder_params = disorder_params if disorder_params is not None else {}
        self.use_ctrw = use_ctrw  # Store flag

        # --- CTRW State ---
        self.wait_times = np.zeros(self.num_walkers, dtype=np.int64)  # Initialize always
        self.ctrw_alpha = self.disorder_params.get('ctrw_alpha', None)
        if self.use_ctrw and (self.ctrw_alpha is None or not (0 < self.ctrw_alpha < 1)):
            raise ValueError("CTRW enabled but 'ctrw_alpha' is missing or invalid.")

        # --- Precomputation ---
        self.precomputed_probs = None
        # Call precomputation if not using default
        if self.disorder_function != self._default_disorder_function:
            print(f"Precomputing probabilities using {self.disorder_function.__name__}...")
            self._precompute_probabilities()
            if self.precomputed_probs is not None:
                print("Precomputation finished.")
            else:
                print("Precomputation failed.")

        self.store_history = store_history
        self.all_positions = []  # Initialize as list
        if self.store_history:
            self.all_positions.append(np.copy(self.initial_positions))


        self.out_of_bounds_walkers = set()

    def _precompute_probabilities(self):
        """
        Calls the vectorized self.disorder_function to precompute probabilities.
        MODIFIED WITH DEBUGGING PRINTS.
        """
        if not callable(self.disorder_function):
            print("Warning: No valid disorder function provided for precomputation.")
            self.precomputed_probs = None
            return

        print(f"Attempting precomputation with {self.disorder_function.__name__}...")
        print(f"  Input xv shape: {self.xv.shape}, dtype: {self.xv.dtype}")
        print(f"  Input yv shape: {self.yv.shape}, dtype: {self.yv.dtype}")
        print(f"  Disorder params: {self.disorder_params}")

        try:
            # --- Step 1: Call the disorder function ---
            print("  Calling disorder function...")
            result_probs = self.disorder_function(
                self.xv, self.yv, **self.disorder_params
            )
            print("  Disorder function call finished.")
            print(
                f"  Function returned type: {type(result_probs)}, shape: {getattr(result_probs, 'shape', 'N/A')}, dtype: {getattr(result_probs, 'dtype', 'N/A')}")

            # --- Step 2: Check for obvious issues before type conversion ---
            if isinstance(result_probs, np.ndarray):
                print(f"  Checking result array for NaNs: {np.isnan(result_probs).any()}")
                print(f"  Checking result array for Infs: {np.isinf(result_probs).any()}")
            else:
                print("  Result is not a NumPy array!")
                raise TypeError("Disorder function did not return a NumPy array.")

            # --- Step 3: Convert type ---
            print(f"  Attempting type conversion to np.float32...")
            self.precomputed_probs = result_probs.astype(np.float32)
            print(
                f"  Type conversion successful. Shape: {self.precomputed_probs.shape}, dtype: {self.precomputed_probs.dtype}")

            # --- Step 4: Check shape ---
            expected_shape = (self.ny, self.nx, 5)
            if self.precomputed_probs.shape != expected_shape:
                raise ValueError(
                    f"Precomputed probs have wrong shape: {self.precomputed_probs.shape}. Expected: {expected_shape}")
            print("  Shape check successful.")

        except Exception as e:
            print(
                f"***** ERROR during precomputation with {self.disorder_function.__name__}: {e} *****")  # Make error stand out
            import traceback
            traceback.print_exc()  # Print the full traceback where the error occurred
            print("Precomputation failed. Check if the disorder function is vectorized correctly.")
            self.precomputed_probs = None  # Ensure it's None if failed

        # Optional sanity check can remain here if desired
        if self.precomputed_probs is not None:
            print("  Running final sum check...")
            sums = np.sum(self.precomputed_probs, axis=2)
            if not np.allclose(sums, 1.0):
                print("  Warning: Precomputed probabilities do not sum to 1 everywhere!")
            else:
                print("  Final sum check passed.")












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


    # --- Modified random_walk_disordered Method ---
    def random_walk_disordered(self):
        """
        Performs one step of disordered walk. Uses standard or CTRW kernel
        based on self.use_ctrw flag. Modifies positions and potentially wait_times.
        """
        if not self.is_uniform_grid:
            raise NotImplementedError("Optimized kernels require uniform grid.")
        if self.precomputed_probs is None:
            raise ValueError("Precomputed probabilities needed but not available.")

        """
        # --- Add Debug Prints ---
        print(f"DEBUG: Inside random_walk_disordered")
        print(f"DEBUG: self.precomputed_probs is None? {self.precomputed_probs is None}")
        print(f"DEBUG: self.disorder_function: {repr(self.disorder_function)}")
        print(f"DEBUG: self._default_disorder_function: {repr(self._default_disorder_function)}")
        print(
            f"DEBUG: self.disorder_function != self._default_disorder_function? {self.disorder_function != self._default_disorder_function}")
        # -----------------------
        """
        if self.use_ctrw:
            # --- Call CTRW Numba kernel ---
            calculated_steps = _run_ctrw_disordered_step_numba(
                self.positions, self.wait_times, self.step, self.precomputed_probs,
                self.x_min_grid, self.y_min_grid, self.dx, self.dy, self.nx, self.ny,
                self.ctrw_alpha
            )
            # wait_times are modified in-place by the kernel
        else:
            # --- Call Standard Numba kernel ---
            calculated_steps = _run_standard_disordered_step_numba(
                self.positions, self.step, self.precomputed_probs,
                self.x_min_grid, self.y_min_grid, self.dx, self.dy, self.nx, self.ny
            )

        # Update positions using the steps returned by the chosen kernel
        self.positions += calculated_steps

        # Apply PBC if needed
        # self.apply_pbc()

    # --------------------------------------------

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

    # In RandomWalk class:
    def random_walk_ordered(self):
        """Calculates one step for the ordered random walk using Numba."""
        steps = random_walk_ordered_numba(self.positions, self.step, self.num_walkers)
        self.positions += steps



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


    # --- trajectories method (ensure wait_times are reset) ---
    def trajectories(self, use_disorder=False):
        self.positions = np.copy(self.initial_positions)
        self.wait_times.fill(0)  # Reset wait times at the start of each trajectory

        self.msd_results = np.zeros(self.num_steps + 1, dtype=np.float64)
        self.msd_results[0] = 0.0

        if self.store_history and not self.all_positions:  # Ensure initial stored if list was cleared
            self.all_positions = [np.copy(self.initial_positions)]


        print(
            f"Running trajectories ({'CTRW' if self.use_ctrw and use_disorder else ('Standard Disordered' if use_disorder else 'Ordered')})...")
        for step_num in range(1, self.num_steps + 1):
            if use_disorder:
                self.random_walk_disordered()  # Will use correct kernel based on self.use_ctrw
            else:
                self.random_walk_ordered()

            if self.store_history:
                self.all_positions.append(np.copy(self.positions))

            displacement = self.positions - self.initial_positions
            sq_displacement = np.sum(displacement ** 2, axis=1)
            self.msd_results[step_num] = np.mean(sq_displacement)


        if self.store_history:
            self.all_positions = np.array(self.all_positions)
        print("Simulation finished. MSD calculated.")
    # ---------------------------------------------------------

    def animate_trajectory(self, walker_index=0, interval=100, save_animation=False, filename="random_walk.gif"):
        """
        Animates the trajectory of a single random walker.
        Handles cases where the walker might not have moved.

        Args:
            walker_index (int): The index of the walker to animate (default: 0).
            interval (int): The delay between frames in milliseconds (default: 100).
            save_animation (bool): Whether to save the animation to a file.
            filename (str): Filename for the saved animation.
        """
        # --- Check if history was stored ---
        if not self.store_history or not hasattr(self, 'all_positions') or len(self.all_positions) == 0:
            print("Error: Cannot animate. Run simulation with store_history=True first.")
            return
        # Ensure all_positions is a numpy array
        if not isinstance(self.all_positions, np.ndarray):
            try:
                # Attempt conversion if it's still a list (should happen at end of trajectories)
                self.all_positions = np.array(self.all_positions)
                print("Converted all_positions list to NumPy array for animation.")
            except Exception as e:
                print(f"Error: Could not convert all_positions to NumPy array: {e}")
                return

        if walker_index >= self.num_walkers:
            print(f"Error: Walker index {walker_index} is out of bounds (0 to {self.num_walkers - 1}).")
            return
        if self.all_positions.shape[0] <= 1:
            print("Error: Not enough time steps in history to animate.")
            return

        fig, ax = plt.subplots()

        # --- Calculate plot limits robustly ---
        x_positions = self.all_positions[:, walker_index, 0]
        y_positions = self.all_positions[:, walker_index, 1]

        x_min, x_max = np.min(x_positions), np.max(x_positions)
        y_min, y_max = np.min(y_positions), np.max(y_positions)

        # Add padding or handle zero range
        x_padding = (x_max - x_min) * 0.1  # 10% padding
        y_padding = (y_max - y_min) * 0.1  # 10% padding

        # If min == max (walker didn't move in that dimension), add default padding
        if np.isclose(x_min, x_max):
            x_padding = 0.1  # Default padding if no movement
            x_min -= x_padding
            x_max += x_padding
        else:
            x_min -= x_padding
            x_max += x_padding

        if np.isclose(y_min, y_max):
            y_padding = 0.1  # Default padding if no movement
            y_min -= y_padding
            y_max += y_padding
        else:
            y_min -= y_padding
            y_max += y_padding

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        # ------------------------------------

        ax.set_xlabel("X Position")
        ax.set_ylabel("Y Position")
        ax.set_title(f"Trajectory of Walker {walker_index}")
        ax.grid(True, which='major', color='grey', linestyle='--')
        # ax.grid(visible=True, which='minor', color='lightgrey', linestyle=':') # Optional minor grid
        # plt.minorticks_on() # Optional minor ticks
        ax.set_aspect('equal', adjustable='box')  # Keep aspect ratio equal

        line, = ax.plot([], [], 'b-', lw=1.5, alpha=0.8)  # Trajectory line
        point, = ax.plot([], [], 'ro', markersize=6)  # Current position marker
        start_point, = ax.plot(x_positions[0], y_positions[0], 'go', markersize=6, label='Start')  # Start marker

        # Ensure update function handles the case of only one frame (though checked above)
        def update(frame):
            if frame == 0:  # Handle first frame explicitly if needed
                line.set_data([], [])
                point.set_data(x_positions[0], y_positions[0])
            else:
                # Plot up to the current frame + 1 (to include the frame itself)
                line.set_data(x_positions[:frame + 1], y_positions[:frame + 1])
                point.set_data(x_positions[frame], y_positions[frame])  # Current point
            return line, point, start_point  # Return all artists being updated

        # Create animation
        # frames should be the number of steps recorded (length of all_positions)
        num_frames = self.all_positions.shape[0]
        ani = animation.FuncAnimation(fig, update, frames=num_frames,
                                      interval=interval, blit=True, repeat=False)

        if save_animation:
            print(f"Attempting to save animation as {filename}...")
            try:
                # Try saving using ffmpeg first (usually better quality)
                ani.save(filename, writer='ffmpeg', fps=1000 / interval, dpi=150)
                print("Animation saved successfully (using ffmpeg).")
            except Exception as e1:
                print(f"  ffmpeg writer failed: {e1}")
                print("  Attempting to save as GIF using pillow...")
                try:
                    # Fallback to pillow for GIF
                    ani.save(filename, writer='pillow', fps=1000 / interval)
                    print("Animation saved successfully (as GIF using pillow).")
                except Exception as e2:
                    print(f"  Pillow writer also failed: {e2}")
                    print("  Could not save animation. Ensure ffmpeg or pillow is installed.")
                    plt.show()  # Show interactively if saving fails
        else:
            plt.show()  # Show interactively if not saving


    def compute_msd(self, direction='all'):
        """
        Compute the mean squared displacement over time.
        MODIFIED: Returns pre-calculated results if available,
                  otherwise calculates from self.all_positions (if stored).
        """
        if hasattr(self, 'msd_results') and direction == 'all':
            # Return the result computed during trajectories() if Option 2 was used
            return self.msd_results
        elif hasattr(self, 'all_positions') and self.all_positions:
            # Fallback to original calculation if all positions were stored
            print("Calculating MSD from stored positions...")
            all_positions_array = np.array(self.all_positions)
            displacement = all_positions_array - self.initial_positions[np.newaxis, :, :]
            # ... (rest of your original MSD calculation for different directions) ...
            if direction == 'x':
                squared_displacement = displacement[:, :, 0] ** 2
            elif direction == 'y':
                squared_displacement = displacement[:, :, 1] ** 2
            elif direction == 'all':
                squared_displacement = np.sum(displacement ** 2, axis=2)
            else:
                raise ValueError("Invalid direction. Choose 'x', 'y', or 'all'.")
            msd = np.mean(squared_displacement, axis=1)
            return msd
        else:
            raise ValueError("Cannot compute MSD. Run trajectories() first, either storing positions or calculating MSD on the fly.")

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
        Dp = self.step ** 2 / (4 * self.dt)  # QUELLO GIUSTO?
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

    def plot_position_histograms(self, time_step, walker_indices=None, num_bins=20):  # Added num_bins
        """
        Plots histograms of walker positions at a given time step and overlays
        a Gaussian fit based on the data's mean and standard deviation.

        Args:
            time_step (int): The time step for which to plot the histograms.
            walker_indices (list, optional): A list of walker indices to include.
                                           If None, all walkers are included.
            num_bins (int): Number of bins to use for the histograms.
        """
        # --- Input Validation ---
        if not self.store_history or not hasattr(self, 'all_positions') or len(self.all_positions) == 0:
            # Updated error message for clarity
            print("Error: Cannot plot histograms. Run simulation with store_history=True first.")
            # Consider raising ValueError instead of just printing
            # raise ValueError("Run the simulation with store_history=True first.")
            return
        if not isinstance(self.all_positions, np.ndarray):
            try:
                self.all_positions = np.array(self.all_positions)
            except Exception as e:
                print(f"Error: Could not convert all_positions to NumPy array: {e}")
                return

        if time_step < 0 or time_step >= self.all_positions.shape[0]:
            print(
                f"Warning: Invalid time step {time_step} for histograms. Max step is {self.all_positions.shape[0] - 1}. Skipping.")
            # Or raise ValueError
            # raise ValueError(f"Invalid time step: {time_step}. Must be between 0 and {self.all_positions.shape[0] - 1}")
            return

        # --- Get Position Data ---
        positions_at_time = self.all_positions[time_step]  # (num_walkers, 2) array
        if walker_indices is None:
            x_positions = positions_at_time[:, 0]
            y_positions = positions_at_time[:, 1]
            num_data_points = self.num_walkers
        else:
            # Ensure indices are valid
            valid_indices = [idx for idx in walker_indices if 0 <= idx < self.num_walkers]
            if not valid_indices:
                print(f"Warning: No valid walker indices provided for histogram at step {time_step}. Skipping.")
                return
            x_positions = positions_at_time[valid_indices, 0]
            y_positions = positions_at_time[valid_indices, 1]
            num_data_points = len(valid_indices)

        if num_data_points == 0:
            print(f"Warning: No data points to plot for histogram at step {time_step}. Skipping.")
            return

        # --- Create Plot ---
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # Slightly wider figure
        fig.suptitle(f'Position Histograms at Time Step {time_step}', fontsize=14)  # Add overall title

        # --- X Histogram and Fit ---
        # Plot histogram and get bin info
        counts_x, bins_x, patches_x = axs[0].hist(x_positions, bins=num_bins, color='skyblue', edgecolor='black',
                                                  alpha=0.7, label='X Data')
        # Calculate statistics for fit
        mean_x = np.mean(x_positions)
        std_x = np.std(x_positions)
        # Create points for the Gaussian curve
        x_fit = np.linspace(bins_x[0], bins_x[-1], 100)
        # Calculate Gaussian PDF
        pdf_x = stats.norm.pdf(x_fit, mean_x, std_x)
        # Scale PDF to match histogram counts (Area under PDF=1, Area under hist=N*bin_width)
        bin_width_x = bins_x[1] - bins_x[0]
        scale_factor_x = num_data_points * bin_width_x
        # Plot scaled Gaussian fit
        axs[0].plot(x_fit, pdf_x * scale_factor_x, 'r--', linewidth=2,
                    label=f'Gaussian Fit\n(μ={mean_x:.2e}, σ={std_x:.2e})')
        axs[0].set_xlabel("X Position")
        axs[0].set_ylabel("Frequency")
        axs[0].set_title("X Positions")
        axs[0].legend()
        axs[0].grid(True, linestyle=':')

        # --- Y Histogram and Fit ---
        # Plot histogram and get bin info
        counts_y, bins_y, patches_y = axs[1].hist(y_positions, bins=num_bins, color='lightgreen', edgecolor='black',
                                                  alpha=0.7, label='Y Data')
        # Calculate statistics for fit
        mean_y = np.mean(y_positions)
        std_y = np.std(y_positions)
        # Create points for the Gaussian curve
        y_fit = np.linspace(bins_y[0], bins_y[-1], 100)
        # Calculate Gaussian PDF
        pdf_y = stats.norm.pdf(y_fit, mean_y, std_y)
        # Scale PDF to match histogram counts
        bin_width_y = bins_y[1] - bins_y[0]
        scale_factor_y = num_data_points * bin_width_y
        # Plot scaled Gaussian fit
        axs[1].plot(y_fit, pdf_y * scale_factor_y, 'r--', linewidth=2,
                    label=f'Gaussian Fit\n(μ={mean_y:.2e}, σ={std_y:.2e})')
        axs[1].set_xlabel("Y Position")
        axs[1].set_ylabel("Frequency")
        axs[1].set_title("Y Positions")
        axs[1].legend()
        axs[1].grid(True, linestyle=':')

        # --- Final Touches ---
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to make room for suptitle
        # plt.show() # Keep this if you want plots to display immediately when called in a loop
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


# =============================================================================
# Multiprocessing setup (Modified run_single_trial)
# =============================================================================
def run_single_trial(params):
    """ Runs one full simulation trial and returns the MSD array. """
    try:
        # Unpack parameters (ensure order matches task_args creation)
        # Renamed disorder_function to disorder_function_param for clarity
        trial_index, num_steps, num_walkers, step, dt, \
        disorder_function_param, disorder_params, xv, yv, use_ctrw_flag, seed = params

        np.random.seed(seed)
        print(f"Starting Trial {trial_index+1} (Seed: {seed}, CTRW: {use_ctrw_flag})...")

        # Instantiate RandomWalk, passing the original disorder function parameter
        rw = RandomWalk(
            num_steps=num_steps, num_walkers=num_walkers, step=step, dt=dt,
            xv=xv, yv=yv, disorder_function=disorder_function_param, # Pass the param here
            disorder_params=disorder_params, use_ctrw=use_ctrw_flag
        )

        # *** CORRECTED LOGIC ***
        # Determine if the disordered walk should be used based on whether
        # a specific disorder function was provided in the parameters.
        should_use_disorder = (disorder_function_param is not None)

        # Pass the correctly determined flag to the trajectories method
        rw.trajectories(use_disorder=should_use_disorder)

        # Compute (or retrieve) MSD results
        msd_result = rw.compute_msd() # Assumes compute_msd retrieves pre-calculated results

        # print(f"Finished Trial {trial_index+1}.")
        return msd_result
    except Exception as e:
        print(f"!!! Error in Trial {trial_index+1}: {e}")
        import traceback; traceback.print_exc()
        return None


# --- main_parallel function (Modified task_args creation) ---
def main_parallel(num_trials_total, num_steps, num_walkers, step, dt,
                  disorder_function, disorder_params, xv, yv, use_ctrw_flag): # Added use_ctrw_flag
    # ... (timer start, get num_workers) ...
    start_time = timer.time(); num_workers = os.cpu_count(); print(f"Detected {num_workers} cores.")

    base_seed = np.random.randint(10000)
    task_args = []
    for i in range(num_trials_total):
        unique_seed = base_seed + i
        task_args.append( # Ensure all needed args are included in the correct order
            (i, num_steps, num_walkers, step, dt,
             disorder_function, disorder_params, xv, yv, use_ctrw_flag, unique_seed) # Added flag
        )

    print(f"\nStarting {num_trials_total} trials using {num_workers} worker processes (CTRW Mode: {use_ctrw_flag})...")
    pool = multiprocessing.Pool(processes=num_workers)
    results = []
    try:
        results = pool.map(run_single_trial, task_args)
    except Exception as e: print(f"!!! Error during parallel execution: {e}")
    finally: pool.close(); pool.join()
    print(f"\nParallel execution finished. Time taken: {timer.time() - start_time:.2f} seconds")

    # ... (Process results as before: filter None, stack, mean) ...
    successful_results = [res for res in results if res is not None]
    if not successful_results: return None, None
    print(f"Successful trials: {len(successful_results)}/{num_trials_total}")
    msd_stack = np.stack(successful_results, axis=0)
    avg_msd = np.mean(msd_stack, axis=0)
    time_axis = np.arange(num_steps + 1)*dt
    return avg_msd, time_axis


# =============================================================================
#IMPORT VECTORIZED DISORDER FUNCTIONS FROM ANOTHER FILE

from vectorized_disorder_funcs import (
    my_spatial_disorder_vectorized,
    gaussian_rest_prob_vectorized,
    uniform_rest_prob_vectorized,
    plateau_rest_prob_vectorized,
    multi_center_rest_prob_vectorized,
    exponential_rest_prob_vectorized,
    boundary_dependent_rest_prob_vectorized
)
# Use the VECTORIZED versions suitable for precomputation
AVAILABLE_DISORDER_FUNCTIONS = {
    "none": None, # Special case for ordered walk
    "uniform": uniform_rest_prob_vectorized,
    "gaussian": gaussian_rest_prob_vectorized,
    "plateau": plateau_rest_prob_vectorized,
    "multi_center": multi_center_rest_prob_vectorized,
    "exponential": exponential_rest_prob_vectorized,
    "boundary": boundary_dependent_rest_prob_vectorized,
    "fixed_rest": my_spatial_disorder_vectorized, # Example name for the fixed rest one
}
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Random Walk Simulation")

    # --- Core Simulation Parameters ---
    parser.add_argument('--steps', type=int, default=1000, help='Number of simulation steps')
    parser.add_argument('--walkers', type=int, default=100, help='Number of walkers')
    parser.add_argument('--trials', type=int, default=10, help='Number of parallel trials')
    parser.add_argument('--step_size', type=float, default=0.001, help='Step size per move')
    parser.add_argument('--dt', type=float, default=0.0001, help='Time step duration')
    parser.add_argument('--grid_size', type=int, default=2000, help='Grid resolution (grid_size x grid_size)')

    # --- CTRW Parameters ---
    parser.add_argument('--ctrw', action='store_true', help='Enable Continuous Time Random Walk (CTRW)')
    parser.add_argument('--alpha', type=float, default=0.7, help='CTRW exponent alpha (0 < alpha < 1)')

    # --- Disorder Function Selection ---
    parser.add_argument('--disorder', type=str, default='none',
                        choices=AVAILABLE_DISORDER_FUNCTIONS.keys(),
                        help='Type of disorder function to use')

    # --- Disorder Function Parameters (add arguments for parameters of ALL functions) ---
    # Note: Only parameters relevant to the chosen --disorder will be used.
    # Gaussian / Exponential / Plateau / Multi-Center / Fixed
    parser.add_argument('--max_rest', type=float,
                        help='Max resting probability (for gaussian, plateau, exp, multi_center, fixed_rest)')
    parser.add_argument('--sigma', type=float, help='Sigma for Gaussian disorder')
    # Uniform
    parser.add_argument('--rest_level', type=float, help='Uniform rest level')
    # Plateau
    parser.add_argument('--plateau_radius', type=float, help='Radius for plateau disorder')
    parser.add_argument('--decay_rate', type=float, help='Decay rate (for plateau, exp, multi_center, boundary)')
    # Multi-Center (simplified example, could add center coords too)
    parser.add_argument('--strength1', type=float, help='Strength for multi-center 1')
    parser.add_argument('--strength2', type=float, help='Strength for multi-center 2')
    # Boundary
    parser.add_argument('--boundary_strength', type=float, help='Strength for boundary disorder')
    # Fixed Rest (my_spatial_disorder_vectorized)
    parser.add_argument('--rest_fixed', type=float, help='Fixed rest probability for fixed_rest type')

    # --- Animation Control ---
    parser.add_argument('--animate', action='store_true', help='Run a single trial and generate animation')
    parser.add_argument('--anim_steps', type=int, default=500, help='Number of steps for animation run')
    parser.add_argument('--anim_walker', type=int, default=0, help='Index of walker to animate')
    parser.add_argument('--save_anim', action='store_true',help='Save the animation file (requires --animate)')
    parser.add_argument('--anim_file', type=str, default='walk_animation.gif', help='Output filename for animation')

    # --- Histogram Control ---
    parser.add_argument('--histograms', action='store_true',
                        help='Run a single trial and show position histograms')
    parser.add_argument('--hist_steps', type=int, nargs='+',  # Expect one or more integers
                        help='List of time steps (integers) to plot histograms for (requires --histograms)')

    args = parser.parse_args()

    # --- Validate Histogram Arguments ---
    if args.histograms and not args.hist_steps:
        parser.error("--histograms requires --hist_steps to be specified.")
    if args.hist_steps and not args.histograms:
        print("Warning: --hist_steps provided but --histograms flag is missing. Histograms will not be generated.")
        # Or parser.error if you want it to be strict

    args = parser.parse_args()

    # --- Select the Disorder Function ---
    selected_disorder_func = AVAILABLE_DISORDER_FUNCTIONS[args.disorder]

    # --- Build Disorder Parameters Dictionary ---
    # Include only non-None arguments relevant to the selected function (or CTRW)
    disorder_params = {}
    potential_params = {
        'max_rest_strength': args.max_rest,  # Name used in gaussian_rest_prob_vectorized
        'sigma': args.sigma,
        'rest_level': args.rest_level,
        'plateau_radius': args.plateau_radius,
        'max_rest': args.max_rest,  # Name used in plateau, exp, multi_center
        'decay_rate': args.decay_rate,
        'strength1': args.strength1,
        'strength2': args.strength2,
        'max_total_rest': args.max_rest,  # Used in multi_center, boundary
        'boundary_strength': args.boundary_strength,
        'rest_fixed': args.rest_fixed,  # Used in my_spatial_disorder_vectorized
    }
    for key, value in potential_params.items():
        if value is not None:
            disorder_params[key] = value

    # Add CTRW alpha if enabled
    if args.ctrw:
        if not (0 < args.alpha < 1):
            parser.error("--alpha must be between 0 and 1 for CTRW")
        disorder_params['ctrw_alpha'] = args.alpha
        print(f"CTRW Enabled with alpha = {args.alpha}")
    else:
        print("CTRW Disabled (Standard Rest/Movement)")

    # --- Setup Grid ---
    print(f"Setting up grid ({args.grid_size}x{args.grid_size})...")
    XV, YV = np.meshgrid(np.linspace(-1, 1, args.grid_size), np.linspace(-1, 1, args.grid_size))
    print("Grid setup done.")

    # --- Run Parallel Simulation for MSD ---
    print(f"\n--- Running Parallel Simulation ({args.trials} Trials) ---")
    print(f"Disorder Function: {args.disorder}")
    print(f"Parameters: {disorder_params}")

    avg_msd, time_axis = main_parallel(
        num_trials_total=args.trials,
        num_steps=args.steps,
        num_walkers=args.walkers,
        step=args.step_size,
        dt=args.dt,
        disorder_function=selected_disorder_func,  # Pass the selected function object
        disorder_params=disorder_params,  # Pass the constructed params
        xv=XV,
        yv=YV,
        use_ctrw_flag=args.ctrw
    )

    # --- Quantitative Analysis ---
    # --- Quantitative Analysis ---
    if avg_msd is not None and time_axis is not None:
        print("\n" + "=" * 30)
        print(" Quantitative Analysis Results")
        print("=" * 30)

        # Define fit range (e.g., last half of the data, avoiding first few points)
        min_fit_step = max(10, args.steps // 2)  # Start fit from step 10 or halfway, whichever is later
        max_fit_step = args.steps
        print(f"Analysis Range Steps: [{min_fit_step}, {max_fit_step}]")  # Print range once

        # Initialize fit results to NaN
        fitted_alpha_msd = np.nan
        fit_intercept_msd = np.nan
        fitted_alpha_msd_t = np.nan

        # Check if the range is valid
        if max_fit_step > min_fit_step and len(time_axis) > max_fit_step:
            # Get indices corresponding to the step range
            idx_min = min_fit_step
            idx_max = max_fit_step
            time_fit_range = time_axis[idx_min: idx_max + 1]
            msd_fit_range = avg_msd[idx_min: idx_max + 1]

            # --- 1. Fit Log-Log MSD to find alpha exponent ---
            print("\n--- Method 1: Log-Log MSD Fit (log(MSD) vs log(t)) ---")
            valid_fit_indices_msd = (time_fit_range > 1e-15) & (msd_fit_range > 1e-15)
            if np.sum(valid_fit_indices_msd) >= 2:
                log_time_msd = np.log(time_fit_range[valid_fit_indices_msd])
                log_msd = np.log(msd_fit_range[valid_fit_indices_msd])
                try:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(log_time_msd, log_msd)
                    fitted_alpha_msd = slope
                    fit_intercept_msd = intercept  # log(C)
                    print(f"  Estimated Alpha (Slope) = {fitted_alpha_msd:.4f}")
                    print(f"  Standard Error          = {std_err:.4f}")
                    print(f"  R-squared               = {r_value ** 2:.4f}")
                except Exception as e:
                    print(f"  Error during log-log MSD fit: {e}")
            else:
                print(f"  Not enough valid data points ({np.sum(valid_fit_indices_msd)}) for log-log MSD fit.")

            # --- 2. Fit Log-Log MSD/Time to find alpha-1 exponent ---
            print("\n--- Method 2: Log-Log MSD/Time Fit (log(MSD/t) vs log(t)) ---")
            msd_over_time_fit_range = msd_fit_range / time_fit_range  # Calculate MSD/t for the range
            valid_fit_indices_msd_t = (time_fit_range > 1e-15) & (msd_over_time_fit_range > 1e-15)  # Check MSD/t > 0
            if np.sum(valid_fit_indices_msd_t) >= 2:
                log_time_msd_t = np.log(time_fit_range[valid_fit_indices_msd_t])
                log_msd_over_time = np.log(msd_over_time_fit_range[valid_fit_indices_msd_t])
                try:
                    slope_alpha_minus_1, intercept_b, r_value_b, p_value_b, std_err_b = stats.linregress(log_time_msd_t,
                                                                                                         log_msd_over_time)
                    # Implied alpha = slope + 1
                    fitted_alpha_msd_t = slope_alpha_minus_1 + 1.0
                    print(f"  Estimated Alpha-1 (Slope) = {slope_alpha_minus_1:.4f}")
                    print(f"  Implied Alpha             = {fitted_alpha_msd_t:.4f}")
                    print(f"  Standard Error (Slope)  = {std_err_b:.4f}")
                    print(f"  R-squared                 = {r_value_b ** 2:.4f}")
                except Exception as e:
                    print(f"  Error during log-log MSD/t fit: {e}")
            else:
                print(f"  Not enough valid data points ({np.sum(valid_fit_indices_msd_t)}) for log-log MSD/t fit.")

        else:
            print(f"  Fit range [{min_fit_step}, {max_fit_step}] invalid or insufficient data length.")



        # --- 3. Calculate Effective Diffusion Coefficient ---
        print("\n--- Effective Diffusion Coefficient (D_eff = MSD / 4t) ---")
        # Calculate D_eff over the same fit range used for alpha
        if max_fit_step > min_fit_step and len(time_axis) > max_fit_step:
            time_eff = time_axis[idx_min: idx_max + 1]
            msd_eff = avg_msd[idx_min: idx_max + 1]
            valid_eff_indices = (time_eff > 1e-15)  # Avoid division by zero

            if np.any(valid_eff_indices):
                # Calculate D_eff = MSD / (4*t) for valid points in the range
                d_eff_values = msd_eff[valid_eff_indices] / (4 * time_eff[valid_eff_indices])

                # Report the average D_eff in the fit range
                avg_d_eff = np.mean(d_eff_values)
                # Report D_eff at the end of the fit range
                final_d_eff = d_eff_values[-1]

                print(f"Analysis Range Steps: [{min_fit_step}, {max_fit_step}]")
                print(f"  Average D_eff in range = {avg_d_eff:.4e}")
                print(f"  Final D_eff in range   = {final_d_eff:.4e}")
            else:
                print("  No valid time points > 0 in range for D_eff calculation.")
                avg_d_eff = np.nan
                final_d_eff = np.nan
        else:
            print(f"  Fit range [{min_fit_step}, {max_fit_step}] invalid or insufficient data.")
            avg_d_eff = np.nan
            final_d_eff = np.nan

        # --- 4. Theoretical Comparison (for ordered case) ---
        # Calculate theoretical D for the standard ordered walk
        D_theory_ordered = args.step_size ** 2 / (4 * args.dt)
        print("\n--- Theoretical Comparison ---")
        print(f"  Theoretical D (Ordered Walk) = {D_theory_ordered:.4e}")
        if args.disorder == 'none' and not args.ctrw:
            print(f"  (Simulation matches theoretical D if Avg D_eff -> Theoretical D)")
        else:
            print(f"  (Disorder/CTRW expected to reduce D_eff compared to theoretical)")

        print("=" * 30 + "\n")

        # --- Plotting ---
        print("--- Plotting Averaged Results ---")
        # Plot MSD/Time
        plt.figure(figsize=(10, 6))
        valid_div = time_axis > 1e-15
        if np.any(valid_div):
            msd_over_time = avg_msd[valid_div] / time_axis[valid_div]
            plt.plot(time_axis[valid_div], msd_over_time, label=f'MSD/Time ({args.disorder})')
            # Add horizontal line for theoretical D*4 (only makes sense for normal diffusion)
            plt.axhline(4 * D_theory_ordered, color='r', linestyle='--', alpha=0.7,
                        label=f'4 * D_theory (Ordered) = {4 * D_theory_ordered:.2e}')
        plt.xlabel('Time (s)')
        plt.ylabel('MSD / Time')
        plt.title(f'Avg Effective Diffusion Coefficient ({args.trials} Trials)')
        plt.grid(True);
        plt.legend();
        plt.show()

        # Plot Log-Log MSD
        plt.figure(figsize=(10, 6))
        valid_log = (time_axis > 1e-15) & (avg_msd > 1e-15)
        if np.any(valid_log):
            plt.loglog(time_axis[valid_log], avg_msd[valid_log], label=f'MSD ({args.disorder})')
            # Plot theoretical line
            slope_1_line = 4 * D_theory_ordered * time_axis[valid_log]
            plt.loglog(time_axis[valid_log], slope_1_line, 'r--', alpha=0.7,
                       label=f'Slope=1 (Theory D={D_theory_ordered:.2e})')
            # Plot the fitted line if fit was successful
            if not np.isnan(fitted_alpha_msd) and not np.isnan(fit_intercept_msd):
                # Calculate fitted line: MSD = exp(intercept) * t^alpha
                fit_line_msd = np.exp(fit_intercept_msd) * (time_axis[valid_log] ** fitted_alpha_msd)
                plt.loglog(time_axis[valid_log], fit_line_msd, 'g:', alpha=0.9, linewidth=2,
                           label=f'Fit (alpha={fitted_alpha_msd:.3f})')

        plt.xlabel('Time (s)')
        plt.ylabel('MSD')
        plt.title(f'Avg Mean Squared Displacement (Log-Log, {args.trials} Trials)')
        plt.grid(True, which='both');
        plt.legend();
        plt.show()

    else:
        print("Parallel simulation failed or produced no results, skipping analysis and plotting.")

        # --- Optional: Run Single Trial for Animation ---
    if args.animate:
        print("\n--- Running Single Trial for Animation ---")
        # Setup grid for animation (can be same or different)
        XV_anim, YV_anim = np.meshgrid(np.linspace(-1, 1, args.grid_size),np.linspace(-1, 1, args.grid_size))  # Use grid_size for consistency

        # Instantiate RandomWalk with store_history=True
        # Make sure RandomWalk class definition exists above
        try:
            rw_anim = RandomWalk(
                num_steps=args.anim_steps,
                num_walkers=args.walkers,  # Use same number of walkers
                step=args.step_size,
                dt=args.dt,
                xv=XV_anim,
                yv=YV_anim,
                disorder_function=selected_disorder_func,  # Use same selected function
                disorder_params=disorder_params,  # Use same constructed params
                use_ctrw=args.ctrw,
                store_history=True  # <<< Enable history storage
            )
        except NameError:
            print("ERROR: RandomWalk class not defined before animation block.")
            # Handle error appropriately, maybe exit
            import sys

            sys.exit(1)
        except Exception as e:
            print(f"ERROR: Failed to initialize RandomWalk for animation: {e}")
            import sys

            sys.exit(1)

        # Determine if disorder is used for this specific run
        use_disorder_anim = (selected_disorder_func is not None)

        # Run the trajectories method
        print(f"Running animation trajectory ({args.anim_steps} steps)...")
        rw_anim.trajectories(use_disorder=use_disorder_anim)

        # Generate animation
        print("Generating/Showing animation...")
        # Make sure animate_trajectory method exists in RandomWalk class
        try:
            rw_anim.animate_trajectory(
                walker_index=args.anim_walker,
                interval=50,  # Example interval
                save_animation=args.save_anim,  # <<< Use the flag value here
                filename=args.anim_file
            )
            if args.save_anim:  # Optional: Print message only if saving attempt was made
                # Note: animate_trajectory should print success/failure messages
                print(f"Animation saving process initiated for {args.anim_file}")
            else:
                print("Animation displayed interactively.")

        except AttributeError:
            print("ERROR: animate_trajectory method not found in RandomWalk class.")
        except Exception as e:
            print(f"ERROR: Failed during animation generation/display: {e}")

        # --- Optional: Run Single Trial for Histograms ---
    if args.histograms and args.hist_steps:
        print("\n--- Running Single Trial for Histograms ---")
        # Determine the maximum step needed for histograms
        max_hist_step = max(args.hist_steps)
        # Ensure the simulation runs long enough for the latest histogram
        hist_run_steps = max(args.steps, max_hist_step)  # Use main steps or max hist step, whichever is longer
        print(f"Running simulation up to step {hist_run_steps} to generate requested histograms.")

        # Setup grid for histogram run (can be same or different)
        XV_hist, YV_hist = np.meshgrid(np.linspace(-1, 1, args.grid_size), np.linspace(-1, 1, args.grid_size))

        # Instantiate RandomWalk with store_history=True
        try:
            rw_hist = RandomWalk(
                num_steps=hist_run_steps,  # Run enough steps
                num_walkers=args.walkers,
                step=args.step_size,
                dt=args.dt,
                xv=XV_hist,
                yv=YV_hist,
                disorder_function=selected_disorder_func,
                disorder_params=disorder_params,
                use_ctrw=args.ctrw,
                store_history=True  # <<< MUST store history
            )
        except NameError:
            print("ERROR: RandomWalk class not defined before histogram block.")
            import sys;

            sys.exit(1)
        except Exception as e:
            print(f"ERROR: Failed to initialize RandomWalk for histograms: {e}")
            import sys;

            sys.exit(1)

        # Determine if disorder is used
        use_disorder_hist = (selected_disorder_func is not None)

        # Run the trajectories method
        print(f"Running histogram trajectory ({hist_run_steps} steps)...")
        rw_hist.trajectories(use_disorder=use_disorder_hist)

        # Generate histograms
        print("Generating Histograms...")
        # Make sure plot_position_histograms method exists
        try:
            for step_to_plot in args.hist_steps:
                if step_to_plot <= hist_run_steps:
                    print(f"  Plotting histogram for step {step_to_plot}...")
                    # Ensure the method exists and handles potential errors
                    rw_hist.plot_position_histograms(time_step=step_to_plot)
                else:
                    # This case shouldn't happen due to hist_run_steps calculation, but good practice
                    print(
                        f"  Warning: Requested histogram step {step_to_plot} exceeds simulation length ({hist_run_steps}). Skipping.")
            # Ensure plots are displayed if running interactively
            plt.show()  # Add this if plots don't show automatically
        except AttributeError:
            print("ERROR: plot_position_histograms method not found in RandomWalk class.")
        except Exception as e:
            print(f"ERROR: Failed during histogram generation: {e}")

    print("\nSimulation Finished.")


