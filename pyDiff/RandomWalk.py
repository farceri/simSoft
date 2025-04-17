import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numba
'''
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
                
UPDATE 16/04: Inserito il resting time. NON TORNA ANCORA IL COEFFIENTE DI DIFFUSIONE CON QUELLO TEORICO

TODO:AGGIUNGERE UN CONTROLLO SULLA POSIZIONE DEI WALKER FUORI DALLA GRIGLIA.
    vedi RESTING-TIME-CLASS-NUMBA-CONTROLBOUNDS.TXT, AGGIUNGERE LE PERIODIC BOUNDARY CONDITIONS. CAPIRE CON GIT
    CAPIRE COME USARE IL PROGRAMMA PER FARE PREDIZIONI FISICHE SUL COEFFICIENTE DI DIFFUSIONE ETC...


'''

'''

GRIGLIA 2D con 2P+2Q+R=1

cambio queste probabilità a seconda del disordine che voglio simulare


'''
class RandomWalk:
    def __init__(self,disorder_function=None,num_steps=1000,step=0.001,num_walkers=100,xv=np.meshgrid(np.linspace(-1,1,2000),np.linspace(-1,1,2000))[0],yv=np.meshgrid(np.linspace(-1,1,2000),np.linspace(-1,1,2000))[1],dt=0.001):
        self.num_walkers = num_walkers
        self.num_steps = num_steps
        self.xv=xv
        self.yv=yv
        self.dt=dt
        self.positions=np.zeros((self.num_walkers,2))
        self.initial_positions=np.copy(self.positions)
        #self.initial_positions=np.zeros((self.num_walkers,2))+0.5
        self.step=step
        self.all_positions = [np.copy(self.positions)]
        self.time = np.arange(self.num_steps + 1)
        self.disorder_function = disorder_function if disorder_function is not None else self._default_disorder_function

    """
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
            #print(f"    probabilities: {probabilities}")

            # Ensure probabilities are valid (basic check)
            if len(probabilities) != 5 or not np.isclose(np.sum(probabilities), 1.0) or np.any(probabilities < 0):
                probabilities = self._default_disorder_function(grid_value_x, grid_value_y)
                #print(f"Warning: Invalid probabilities at ({grid_value_x:.3f}, {grid_value_y:.3f}), using default.")

            rest_prob = probabilities[4]
            # print(f"    rest_prob: {rest_prob}")# Get resting probability from function


           
            rand_val=np.random.rand()
            #print(f"    rand_val: {rand_val}")
            if rand_val >= rest_prob:
                #print("    Walker is moving")
                #moved = True
                direction_choice = np.random.choice(4, p=probabilities[0:4])
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
        """
        Simulates a disordered random walk for all walkers in a single time step.
        """
        steps = np.zeros((self.num_walkers, 2))

        for i in range(self.num_walkers):
            y_idx, x_idx = self._get_grid_index(self.positions[i, 0], self.positions[i, 1])
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
                    #direction_choice = np.random.choice(4, p=probabilities[:4])

                if direction_choice == 0:  # +x
                    steps[i, 0] = self.step
                elif direction_choice == 1:  # -x
                    steps[i, 0] = -self.step
                elif direction_choice == 2:  # +y
                    steps[i, 1] = self.step
                elif direction_choice == 3:  # -y
                    steps[i, 1] = -self.step

                self.positions[i] += steps[i]  # Update position for THIS walker ONLY

        return self.positions



    def _default_disorder_function(self, x, y):
        """Default: Uniform probability distribution (no spatial disorder)."""
        return np.array([0.25, 0.25, 0.25, 0.25,0]) # [+x, -x, +y, -y, rest]

    def _get_grid_index(self, x, y):
        """Find the indices of the nearest grid point."""
        x_index = np.argmin(np.abs(np.linspace(-1, 1, self.xv.shape[1]) - x))
        y_index = np.argmin(np.abs(np.linspace(-1, 1, self.yv.shape[0]) - y))
        return y_index, x_index



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


    def trajectories(self, use_disorder=False):
        self.all_positions = [np.copy(self.initial_positions)]  # Reset positions
        for _ in range(self.num_steps):
            if use_disorder:
                self.random_walk_disordered()
            else:
                self.random_walk_ordered()
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

    def animate_trajectory(self, walker_index=0, interval=100,save_animation=False, filename="random_walk.gif"):
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
        point, = ax.plot([], [], 'ro', markersize=8) # Current position

        def update(frame):
            xdata = self.all_positions[:frame, walker_index, 0]
            ydata = self.all_positions[:frame, walker_index, 1]
            line.set_data(xdata, ydata)
            point.set_data(xdata[-1:], ydata[-1:]) # Update current position marker
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
        Dp=self.step**2/(2*self.dt) #QUELLO GIUSTO?
        D=0.25*msd_valid/time_valid

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
        return Dp,D,slope,intercept


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
    rest_prob =0.1 + 0.2 * np.abs(x * y)  # Example: Higher rest probability near origin

    '''
    #
    '''
    rest_prob = 0.5
    prob_plus_x = 0.25*(1-rest_prob)
    prob_minus_x = 0.25*(1-rest_prob)
    prob_plus_y = 0.25*(1-rest_prob)
    prob_minus_y =0.25*(1-rest_prob)
    '''
    prob_plus_x = 1 / np.sqrt(np.pi * (2 * 0.1)) * np.exp(-x ** 2 / (2 * 0.1))
    prob_minus_x = 1 / np.sqrt(np.pi * (2 * 0.1)) * np.exp(-x ** 2 / (2 * 0.1))
    prob_plus_y = 1 / np.sqrt(np.pi * (2 * 0.1)) * np.exp(-y ** 2 / (2 * 0.1))
    prob_minus_y = 1 / np.sqrt(np.pi * (2 * 0.1)) * np.exp(-y ** 2 / (2 * 0.1))
    rest_prob = 1 / np.sqrt(np.pi * (4 * 0.5)) * np.exp(-(y ** 2+x**2) / (4 * 0.5))






    probs = np.clip([prob_plus_x, prob_minus_x, prob_plus_y, prob_minus_y,rest_prob], 0, 1)
    #probs=[prob_plus_x, prob_minus_x, prob_plus_y, prob_minus_y,rest_prob]

    if np.sum(probs) > 0:
        probs = probs / np.sum(probs)
    else:
        probs = [0.25, 0.25, 0.25, 0.25,0]

    return probs


def main():
    rw = RandomWalk(disorder_function=my_spatial_disorder)
    #rw=RandomWalk()
    time = np.arange(rw.num_steps + 1)


    rw.trajectories(use_disorder=True)
    #rw.trajectories(use_disorder=False)



    #rw.plot_position_histograms(time_step=100)



    total_msd=rw.compute_msd()
    x_msd=rw.compute_msd(direction='x')
    y_msd=rw.compute_msd(direction='y')


    D=rw.compute_D()
    print("Dp=",D[0])
    print("slope (no disorder)[2]:",D[2])
    print("intercept (no disorder)[]:", np.exp(D[3])/4)
    

    # Plot the trajectory of the first walker
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

    '''
    fig, ax = plt.subplots()
    valid_indices = (time > 0)
    time_valid = time[valid_indices]
    ax.plot(time_valid,D[1])
    plt.show()
    '''
    #rw.animate_trajectory(walker_index=0, interval=50)



if __name__ == "__main__":
    main()