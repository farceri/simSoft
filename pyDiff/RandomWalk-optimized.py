import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numba
import multiprocessing # Import the module
import time as timer # To time the execution
import os # To potentially get CPU count
import argparse # Import argparse
import yaml
import traceback
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
    """ Draws wait time from tau^-(1+alpha). alpha MUST be > 0. """
    # Add check for safety, although clipping should prevent alpha=0
    if alpha <= 0:
         # Return minimum wait time if alpha is invalid
         # Or handle differently (e.g., raise error if Numba supported it well here)
         return np.int64(1)
    v = 1.0 - np.random.rand() # Random number in (0, 1]
    # Calculate wait steps: (1-v)^(-1/alpha) is equivalent for uniform v in [0,1)
    # Using ceil ensures wait_steps >= 1
    wait_steps = np.int64(np.ceil(v**(-1.0 / alpha)))
    # Numba might warn about potential division by zero if alpha could be 0
    # Clipping alpha beforehand is important.
    return max(np.int64(1), wait_steps) # Ensure minimum wait is 1 step

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
    precomputed_probs,   # Input: The table of probabilities for each grid site (ny, nx, 5)
    alpha_grid,          # <<< INPUT: Grid of alpha values (ny, nx) >>>
    x_min_grid, y_min_grid, dx, dy, nx, ny # Input: Grid info for indexing
    ):
    """ Numba kernel for one CTRW step with position-dependent power-law waits. """
    num_walkers = positions.shape[0]
    steps = np.zeros_like(positions)

    for i in range(num_walkers):
        if wait_times[i] > 0:
            wait_times[i] -= 1
        else:
            current_x = positions[i, 0]
            current_y = positions[i, 1]
            y_idx, x_idx = _calculate_grid_index_fast_numba(
                current_x, current_y, x_min_grid, y_min_grid, dx, dy, nx, ny)

            probabilities = precomputed_probs[y_idx, x_idx, :]
            rest_prob = probabilities[4]

            if np.random.rand() < rest_prob:
                # --- Walker CHOOSES TO REST ---
                # <<< Look up local alpha value from the grid >>>
                local_alpha = alpha_grid[y_idx, x_idx]

                # <<< Draw wait time using the local alpha >>>
                tau = draw_power_law_wait_time(local_alpha)

                wait_times[i] = max(np.int64(0), tau - 1)
            else:
                # --- Walker CHOOSES TO MOVE ---
                direction_probs = probabilities[:4]
                direction_probs_sum = np.sum(direction_probs)
                if direction_probs_sum > 1e-9:
                    choice_rand = np.random.rand() * direction_probs_sum
                    p_plus_x  = direction_probs[0]
                    p_minus_x = direction_probs[1]
                    p_plus_y  = direction_probs[2]

                    if choice_rand < p_plus_x: direction_choice = 0
                    elif choice_rand < (p_plus_x + p_minus_x): direction_choice = 1
                    elif choice_rand < (p_plus_x + p_minus_x + p_plus_y): direction_choice = 2
                    else: direction_choice = 3

                    if direction_choice == 0: steps[i, 0] = step_size
                    elif direction_choice == 1: steps[i, 0] = -step_size
                    elif direction_choice == 2: steps[i, 1] = step_size
                    elif direction_choice == 3: steps[i, 1] = -step_size
                # else: Cannot move if rest_prob is 1.0

    return steps

@numba.njit(cache=True)
def _run_direction_disordered_step_numba(
    positions,           # Input: Current (x,y) for all walkers
    step_size,           # Input: How far a walker moves
    precomputed_probs,   # INPUT: The table of probabilities (ny, nx, 4) [P+x, P-x, P+y, P-y]
    x_min_grid, y_min_grid, dx, dy, nx, ny # Input: Grid info for indexing
    ):
    """ Numba kernel for one step with 4 direction probabilities (Prest=0). """
    num_walkers = positions.shape[0]
    steps = np.zeros_like(positions)

    for i in range(num_walkers):
        current_x = positions[i, 0]
        current_y = positions[i, 1]
        y_idx, x_idx = _calculate_grid_index_fast_numba(
            current_x, current_y, x_min_grid, y_min_grid, dx, dy, nx, ny)

        # Probabilities for [P+x, P-x, P+y, P-y]
        direction_probs = precomputed_probs[y_idx, x_idx, :]
        # No resting check needed, sum should be 1.0

        # Choose direction based on the 4 probabilities
        # (Using cumulative sum method - efficient in Numba)
        choice_rand = np.random.rand() # Random number in [0, 1)
        p_plus_x  = direction_probs[0]
        p_minus_x = direction_probs[1]
        p_plus_y  = direction_probs[2]
        # p_minus_y = direction_probs[3] # Not needed for cumulative check

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

# =============================================================================
#IMPORT VECTORIZED DISORDER FUNCTIONS FROM ANOTHER FILE

from vectorized_disorder_funcs import (
    my_spatial_disorder_vectorized,
    gaussian_rest_prob_vectorized,
    uniform_rest_prob_vectorized,
    plateau_rest_prob_vectorized,
    multi_center_rest_prob_vectorized,
    exponential_rest_prob_vectorized,
    boundary_dependent_rest_prob_vectorized,
    constant_alpha,
    gaussian_alpha,
    linear_gradient_alpha,
    biased_towards_origin,
    vortex_flow
)
# Use the VECTORIZED versions suitable for precomputation
AVAILABLE_STANDARD_DISORDER_FUNCTIONS = {
    "none": None, # Special case for ordered walk
    "uniform": uniform_rest_prob_vectorized,
    "gaussian": gaussian_rest_prob_vectorized,
    "plateau": plateau_rest_prob_vectorized,
    "multi_center": multi_center_rest_prob_vectorized,
    "exponential": exponential_rest_prob_vectorized,
    "boundary": boundary_dependent_rest_prob_vectorized,
    "fixed_rest": my_spatial_disorder_vectorized,
    # Example name for the fixed rest one
}

AVAILABLE_ALPHA_FUNCTIONS = {
    "constant": constant_alpha,
    "gaussian": gaussian_alpha,
    "linear_gradient": linear_gradient_alpha,
    # Add keys matching your alpha function names
}

AVAILABLE_DIRECTION_DISORDER_FUNCTIONS = {
    "bias_origin": biased_towards_origin,
    "vortex": vortex_flow,
    # ... other directional functions
}

# Combine for argparse, but keep separate for logic
ALL_AVAILABLE_DISORDER_FUNCTIONS = {
    "none": None,
    **AVAILABLE_STANDARD_DISORDER_FUNCTIONS,
    **AVAILABLE_DIRECTION_DISORDER_FUNCTIONS
}


class RandomWalk:
    def __init__(self,
                 num_steps=1000,
                 num_walkers=100,
                 step=0.001,
                 dt=0.0001,
                 xv=None,  # Provide default grid creation if needed
                 yv=None,
                 disorder_function=None,
                 disorder_params=None,
                 disorder_mode='standard',  # <<< Default value if not passed
                 use_ctrw=False,
                 alpha_type='constant',
                 alpha_params=None,
                 use_pbc=False,
                 check_bounds=False,
                 store_history=False,
                 interpolation=False):  # Added interpolation back if needed

        self.num_walkers = num_walkers
        self.num_steps = num_steps

        # --- Grid Setup ---
        if xv is None or yv is None:
            print("Warning: No grid provided, using default 100x100 grid.")
            xv, yv = np.meshgrid(np.linspace(-1, 1, 100), np.linspace(-1, 1, 100))
        self.xv = xv
        self.yv = yv
        self.ny, self.nx = self.xv.shape
        self.x_min = np.min(self.xv)
        self.x_max = np.max(self.xv)
        self.y_min = np.min(self.yv)
        self.y_max = np.max(self.yv)
        self.box_width = self.x_max - self.x_min
        self.box_height = self.y_max - self.y_min

        # --- Fast Indexing Setup ---
        self.is_uniform_grid = False
        self.x_min_grid, self.y_min_grid, self.dx, self.dy = 0.0, 0.0, 1.0, 1.0
        if self.nx > 1 and self.ny > 1:
            self.x_coords = self.xv[0, :]
            self.y_coords = self.yv[:, 0]
            dxs = np.diff(self.x_coords)
            dys = np.diff(self.y_coords)
            if np.allclose(dxs, dxs[0]) and np.allclose(dys, dys[0]):
                self.is_uniform_grid = True
                self.x_min_grid = self.x_coords[0];
                self.y_min_grid = self.y_coords[0]
                self.dx = dxs[0] if abs(dxs[0]) > 1e-15 else 1.0
                self.dy = dys[0] if abs(dys[0]) > 1e-15 else 1.0
                print("Uniform grid detected. Using fast indexing.")
            else:
                print("Non-uniform grid detected.")

        # --- Simulation Parameters ---
        self.dt = dt
        self.step = step
        self.time = np.arange(self.num_steps + 1) * self.dt
        self.positions = np.zeros((self.num_walkers, 2), dtype=np.float64)
        self.initial_positions = np.copy(self.positions)
        self.store_history = store_history
        self.all_positions = []
        if self.store_history: self.all_positions.append(np.copy(self.initial_positions))

        # --- Boundary Conditions ---
        self.use_pbc = use_pbc
        self.perform_bounds_check = check_bounds
        if self.use_pbc and self.perform_bounds_check:
            print("Warning: Both PBC and bounds checking enabled.")
        self.out_of_bounds_walkers = set()

        # --- Disorder Setup ---
        self.disorder_function = disorder_function
        self.disorder_params = disorder_params if disorder_params is not None else {}
        self.disorder_mode = disorder_mode  # <<< STORE the passed mode (Removed overwrite)
        self.precomputed_probs = None

        # --- CTRW Setup ---
        self.use_ctrw = use_ctrw
        self.wait_times = np.zeros(self.num_walkers, dtype=np.int64)
        self.alpha_params = alpha_params if alpha_params is not None else {}
        self.alpha_grid = None
        self.alpha_function = None
        if self.use_ctrw:
            # Ensure AVAILABLE_ALPHA_FUNCTIONS is accessible here (global or imported)
            self.alpha_function = AVAILABLE_ALPHA_FUNCTIONS.get(alpha_type)
            if not callable(self.alpha_function):
                print(f"ERROR in __init__: Could not find alpha function for type '{alpha_type}'")
                # Decide how to handle: raise error? set self.use_ctrw = False?
        print(f"  [RW Init Debug Trial {os.getpid()}] Received alpha_function: {repr(self.alpha_function)}")
        print(f"  [RW Init Debug Trial {os.getpid()}] Is callable? {callable(self.alpha_function)}")
        # --- Precomputation Logic (Correctly Indented) ---
        # Disable CTRW if mode requires it
        if self.disorder_mode == 'direction_only' and self.use_ctrw:
            print("Warning: CTRW disabled because 'direction_only' disorder mode selected.")
            self.use_ctrw = False
        if self.disorder_mode == 'none' and self.use_ctrw:
            print("Warning: CTRW disabled because disorder mode is 'none'.")
            self.use_ctrw = False

        # Call appropriate probability precomputation based on mode
        if self.disorder_mode == 'standard':
            print("Precomputing standard (5) probabilities...")
            # Ensure this method exists and handles errors
            self._precompute_standard_probabilities()
        elif self.disorder_mode == 'direction_only':
            print("Precomputing directional (4) probabilities...")
            # Ensure this method exists and handles errors
            self._precompute_direction_probabilities()
        elif self.disorder_mode == 'none':
            print("No disorder function provided. Using ordered walk.")
        else:
            # This case should ideally not be reached if mode is validated earlier
            print(f"Warning: Unknown disorder mode '{self.disorder_mode}' during initialization.")
            self.disorder_mode = 'none'  # Fallback to ordered? Or raise error?

        # Precompute alpha grid only if CTRW is still enabled (must be standard mode)
        if self.use_ctrw:  # Check if CTRW is still enabled after potential disabling above
            print("Precomputing alpha grid for CTRW...")
            # Ensure this method exists and handles errors
            self._precompute_alpha_grid()

        # --- Other Attributes ---
        self.msd_results = None
        self.current_step_num = 0  # For error messages



    def _precompute_alpha_grid(self):
        """ Calls the vectorized self.alpha_function to precompute alpha grid. """
        if not callable(self.alpha_function):
            print("Warning: No valid alpha function provided for CTRW alpha grid precomputation.")
            self.alpha_grid = None
            return

        print(f"  Attempting alpha precomputation with {self.alpha_function.__name__}...")
        print(f"  Alpha params: {self.alpha_params}")

        try:
            # Call the assigned alpha function (which must be vectorized)
            result_alpha = self.alpha_function(
                self.xv, self.yv, **self.alpha_params
            )
            # Basic validation
            if not isinstance(result_alpha, np.ndarray):
                raise TypeError("Alpha function did not return a NumPy array.")
            if result_alpha.shape != (self.ny, self.nx):
                raise ValueError(f"Alpha grid shape mismatch: {result_alpha.shape}. Expected {(self.ny, self.nx)}")
            if np.isnan(result_alpha).any() or np.isinf(result_alpha).any():
                print("Warning: NaNs or Infs detected in precomputed alpha grid!")
            # Check if values are reasonable (e.g., mostly between 0 and 1)
            if np.any(result_alpha <= 0) or np.any(result_alpha >= 1):
                print(
                    "Warning: Alpha grid contains values outside typical (0, 1) range after function call (clipping applied).")

            # Ensure correct dtype (e.g., float32 for consistency with Numba)
            self.alpha_grid = result_alpha.astype(np.float32)
            print(
                f"  Alpha grid type conversion successful. Shape: {self.alpha_grid.shape}, dtype: {self.alpha_grid.dtype}")

        except Exception as e:
            print(f"***** ERROR during alpha precomputation: {e} *****")
            import traceback
            traceback.print_exc()
            self.alpha_grid = None  # Ensure it's None if failed







    def _precompute_standard_probabilities(self):
        """ Calls the vectorized self.disorder_function to precompute 5 probabilities. """
        if not callable(self.disorder_function) or self.disorder_function == self._default_disorder_function:
             # Handle default case if needed, maybe precompute uniform [0.25 ... 0]
             print("Using default or no standard disorder function.")
             # If default needed, create the 5-prob array here
             # self.precomputed_probs = self._default_disorder_function(self.xv, self.yv).astype(np.float32)
             return # Or precompute default

        print(f"Attempting standard precomputation with {self.disorder_function.__name__}...")
        try:
            result_probs = self.disorder_function(self.xv, self.yv, **self.disorder_params)
            # --- Validation for 5 probabilities ---
            if not isinstance(result_probs, np.ndarray): raise TypeError("Function didn't return NumPy array.")
            if result_probs.shape != (self.ny, self.nx, 5): raise ValueError(f"Expected shape (..., 5), got {result_probs.shape}")
            # Add NaN/Inf checks
            # Add sum check (should sum to 1.0)
            sums = np.sum(result_probs, axis=-1)
            if not np.allclose(sums, 1.0): print("Warning: Standard probabilities do not sum to 1!")
            self.precomputed_probs = result_probs.astype(np.float32)
            print("Standard precomputation finished.")
        except Exception as e:
            print(f"***** ERROR during standard precomputation: {e} *****")
            self.precomputed_probs = None; raise # Re-raise

    def _precompute_direction_probabilities(self):
        """ Calls the vectorized self.disorder_function to precompute 4 probabilities. """
        if not callable(self.disorder_function):
             print("Error: No valid directional disorder function provided.")
             self.precomputed_probs = None; raise ValueError("Missing directional function")

        print(f"Attempting directional precomputation with {self.disorder_function.__name__}...")
        try:
            result_probs = self.disorder_function(self.xv, self.yv, **self.disorder_params)
            # --- Validation for 4 probabilities ---
            if not isinstance(result_probs, np.ndarray): raise TypeError("Function didn't return NumPy array.")
            if result_probs.shape != (self.ny, self.nx, 4): raise ValueError(f"Expected shape (..., 4), got {result_probs.shape}")
            # Add NaN/Inf checks
            # Add sum check (should sum to 1.0)
            sums = np.sum(result_probs, axis=-1)
            if not np.allclose(sums, 1.0): print("Warning: Directional probabilities do not sum to 1!")
            self.precomputed_probs = result_probs.astype(np.float32)
            print("Directional precomputation finished.")
        except Exception as e:
            print(f"***** ERROR during directional precomputation: {e} *****")
            self.precomputed_probs = None; raise # Re-raise

    def apply_pbc(self):
        """Apply periodic boundary conditions to keep particles inside the simulation box defined by grid min/max."""
        # Wrap X coordinates
        self.positions[:, 0] = self.x_min + (self.positions[:, 0] - self.x_min) % self.box_width
        # Wrap Y coordinates
        self.positions[:, 1] = self.y_min + (self.positions[:, 1] - self.y_min) % self.box_height
        # --- Old version using box_size array (less direct for non-origin centered box) ---
        # self.positions = (self.positions + self.box_size / 2) % self.box_size - self.box_size / 2 # Assumes origin-centered box?
        # --- Alternative if box is centered at 0 and goes from -L/2 to L/2 ---
        # L_x = self.box_width
        # L_y = self.box_height
        # self.positions[:, 0] = (self.positions[:, 0] + L_x / 2) % L_x - L_x / 2
        # self.positions[:, 1] = (self.positions[:, 1] + L_y / 2) % L_y - L_y / 2

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
        """ Performs one step of disordered walk based on the configured mode. """
        if not self.is_uniform_grid:
            raise NotImplementedError("Optimized kernels require uniform grid.")

        # Check if precomputation succeeded for the chosen mode
        if self.precomputed_probs is None:
             # This check might be redundant if precomputation raises errors, but good safety.
             raise ValueError(f"Precomputed probabilities (mode: {self.disorder_mode}) needed but not available.")

        # --- Call the correct kernel based on mode ---
        if self.disorder_mode == 'direction_only':
            calculated_steps = _run_direction_disordered_step_numba(
                self.positions, self.step, self.precomputed_probs,
                self.x_min_grid, self.y_min_grid, self.dx, self.dy, self.nx, self.ny
            )
        elif self.disorder_mode == 'standard':
            if self.use_ctrw:
                if self.alpha_grid is None: raise ValueError("CTRW enabled, but alpha grid missing.")
                calculated_steps = _run_ctrw_disordered_step_numba(
                    self.positions, self.wait_times, self.step, self.precomputed_probs,
                    self.alpha_grid, # Pass alpha grid
                    self.x_min_grid, self.y_min_grid, self.dx, self.dy, self.nx, self.ny
                )
            else:
                calculated_steps = _run_standard_disordered_step_numba(
                    self.positions, self.step, self.precomputed_probs,
                    self.x_min_grid, self.y_min_grid, self.dx, self.dy, self.nx, self.ny
                )
        else: # Should only be 'none' if called incorrectly
             raise ValueError(f"random_walk_disordered called with invalid mode: {self.disorder_mode}")


        # Update positions and apply BCs
        self.positions += calculated_steps
        # if self.use_pbc: self.apply_pbc() ... etc

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
        """ Checks if any walkers are out of the grid and prints a warning. """
        # Use self.x_min, self.x_max etc defined in __init__
        out_of_bounds = np.where(
            (positions[:, 0] < self.x_min) | (positions[:, 0] > self.x_max) |
            (positions[:, 1] < self.y_min) | (positions[:, 1] > self.y_max)
        )[0]
        newly_out = set(out_of_bounds) - self.out_of_bounds_walkers # Find only newly out-of-bounds walkers
        if newly_out: # Only print warning for newly out walkers
            print(f"Warning: Walkers went out of bounds at step {self.current_step_num}: {sorted(list(newly_out))}") # Requires tracking step num
            # print("Out-of-bounds walker indices:", out_of_bounds) # Original print
            self.out_of_bounds_walkers.update(newly_out) # Update the set


    # --- trajectories method (ensure wait_times are reset) ---
    def trajectories(self, use_disorder=False):
        self.positions = np.copy(self.initial_positions)
        self.wait_times.fill(0)  # Reset wait times at the start of each trajectory

        self.msd_results = np.zeros(self.num_steps + 1, dtype=np.float64)
        self.msd_results[0] = 0.0
        self.out_of_bounds_walkers.clear()

        if self.store_history and not self.all_positions:  # Ensure initial stored if list was cleared
            self.all_positions = [np.copy(self.initial_positions)]

        if self.disorder_mode == 'none':
            run_label = "Ordered"
        elif self.disorder_mode == 'direction_only':
            run_label = "Directional Disordered"
        elif self.disorder_mode == 'standard':
            run_label = "Standard Disordered (CTRW)" if self.use_ctrw else "Standard Disordered"
        else:
            run_label = "Unknown"
        print(f"Running trajectories ({run_label})...")


        for step_num in range(1, self.num_steps + 1):
            # --- Choose walk step based on mode ---
            if self.disorder_mode == 'none':
                self.random_walk_ordered()
            elif self.disorder_mode == 'standard' or self.disorder_mode == 'direction_only':
                # Call the unified disordered method, it will use the correct kernel
                self.random_walk_disordered()
            else:
                raise ValueError(f"Invalid disorder mode in trajectories loop: {self.disorder_mode}")

            # --- Apply Boundary Conditions / Checks ---
            if self.use_pbc:
                self.apply_pbc()
            if self.perform_bounds_check:
                self._check_out_of_bounds(self.positions)

            if self.store_history:
                self.all_positions.append(np.copy(self.positions))

            displacement = self.positions - self.initial_positions
            sq_displacement = np.sum(displacement ** 2, axis=1)
            self.msd_results[step_num] = np.mean(sq_displacement)


        if self.store_history:
            self.all_positions = np.array(self.all_positions)
        print("Simulation finished. MSD calculated.")
        if self.perform_bounds_check and self.out_of_bounds_walkers:
            print(f"Total unique walkers that went out of bounds: {len(self.out_of_bounds_walkers)}")

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
    trial_index = -1
    try:
        # --- Unpack parameters including type names and param dicts ---
        trial_index, num_steps, num_walkers, step, dt, \
        disorder_type_param, disorder_params_param, \
        alpha_type_param, alpha_params_param, \
        disorder_mode_param, \
        xv, yv, \
        use_ctrw_flag, use_pbc_flag, check_bounds_flag, seed = params

        np.random.seed(seed)
        print(f"Starting Trial {trial_index+1} (Seed: {seed}, Mode: {disorder_mode_param}, CTRW: {use_ctrw_flag}, DisType: {disorder_type_param}, AlphaType: {alpha_type_param})...")

        # *** Look up the actual function objects based on the type names ***
        # Ensure the dictionaries are accessible in this scope (e.g., global or imported)
        selected_disorder_func_obj = ALL_AVAILABLE_DISORDER_FUNCTIONS.get(disorder_type_param)
        selected_alpha_func_obj = None
        if use_ctrw_flag:
            selected_alpha_func_obj = AVAILABLE_ALPHA_FUNCTIONS.get(alpha_type_param)
            # Add error check if function not found and CTRW is True
            if selected_alpha_func_obj is None and alpha_type_param != 'N/A':
                 raise ValueError(f"Alpha function type '{alpha_type_param}' requested but not found in dictionary.")

        # Instantiate RandomWalk, passing the function OBJECTS
        # Ensure RandomWalk.__init__ accepts these keyword arguments
        rw = RandomWalk(
            num_steps=num_steps, num_walkers=num_walkers, step=step, dt=dt,
            xv=xv, yv=yv,
            disorder_function=selected_disorder_func_obj,  # Pass function object
            disorder_params=disorder_params_param,
            disorder_mode=disorder_mode_param,
            use_ctrw=use_ctrw_flag,
            # === MODIFIED LINES START ===
            alpha_type=alpha_type_param,  # Pass the type string (e.g., 'constant')
            alpha_params=alpha_params_param,  # Pass the parameter dictionary
            # === MODIFIED LINES END ===
            use_pbc=use_pbc_flag,
            check_bounds=check_bounds_flag
            # store_history=...
        )

        # trajectories method now uses the internal disorder_mode
        rw.trajectories()

        msd_result = rw.compute_msd()
        return msd_result
    except Exception as e:
        trial_label = trial_index + 1 if trial_index != -1 else 'UNKNOWN'
        print(f"!!! Error in Trial {trial_label}: {e}")
        traceback.print_exc()
        return None



# --- main_parallel function (Modified task_args creation) ---
def main_parallel(num_trials_total, num_steps, num_walkers, step, dt,
                  disorder_type, disorder_params, # Accept type name & params
                  alpha_type, alpha_params,       # Accept type name & params
                  disorder_mode,
                  xv, yv,
                  use_ctrw_flag, use_pbc_flag, check_bounds_flag):
    start_time = timer.time()
    # ... (get num_workers) ...
    try: num_workers = os.cpu_count(); print(f"Detected {num_workers} cores.")
    except: num_workers=1

    base_seed = np.random.randint(10000)
    task_args = []
    for i in range(num_trials_total):
        unique_seed = base_seed + i
        # Add disorder_mode to the arguments passed to each worker
        # *** CORRECTED: Packing 16 items into the tuple ***
        task_args.append(
            (i, num_steps, num_walkers, step, dt,               # 5
             disorder_type, disorder_params,                    # 2
             alpha_type, alpha_params,                          # 2
             disorder_mode,                                     # 1
             xv, yv,                                            # 2
             use_ctrw_flag, use_pbc_flag, check_bounds_flag,    # 3
             unique_seed)                                       # 1 --> Total 16 items
        )

    # Print the mode being used for the parallel run
    print(f"\nStarting {num_trials_total} trials using {num_workers} worker processes (Disorder Mode: {disorder_mode}, CTRW: {use_ctrw_flag})...")
    # ... (multiprocessing pool execution, result processing) ...
    results = []
    try:
        # Ensure Pool is managed correctly (e.g., with 'with' statement)
        with multiprocessing.Pool(processes=num_workers) as pool:
            results = pool.map(run_single_trial, task_args)
    except Exception as e:
        # This catches errors during the map process itself (like pickling issues)
        print(f"!!! Error during parallel execution setup/map: {e}")
        # Optionally print traceback here too if needed
        # traceback.print_exc()

    print(f"\nParallel execution finished. Time taken: {timer.time() - start_time:.2f} seconds")

    # ... (rest of the function: process results, return avg_msd, time_axis) ...
    successful_results = [res for res in results if res is not None]
    if not successful_results:
        print("No trials completed successfully.")
        return None, None
    print(f"Successful trials: {len(successful_results)}/{num_trials_total}")
    try:
        msd_stack = np.stack(successful_results, axis=0)
        avg_msd = np.mean(msd_stack, axis=0)
        time_axis = np.arange(num_steps + 1) * dt
        return avg_msd, time_axis
    except Exception as e:
        print(f"!!! Error processing results (e.g., stacking): {e}")
        traceback.print_exc()
        return None, None




if __name__ == "__main__":


    # --- Argument Parser: Only for the config file path ---
    parser = argparse.ArgumentParser(description="Run Random Walk Simulation from Config File")
    parser.add_argument('config_file', type=str,
                        help='Path to the YAML configuration file')
    args = parser.parse_args()  # args now only contains args.config_file

    # --- Load Configuration from YAML File ---
    try:
        with open(args.config_file, 'r') as f:
            config = yaml.safe_load(f)
        print(f"Configuration loaded successfully from {args.config_file}")
    except FileNotFoundError:
        print(f"Error: Configuration file not found at {args.config_file}")
        import sys;

        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file {args.config_file}: {e}")
        import sys;

        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred loading config: {e}")
        import sys;

        sys.exit(1)

    # --- Extract Parameters (using .get() for safety/defaults) ---
    # Simulation Params
    sim_params = config.get('simulation', {})
    steps = sim_params.get('steps', 1000)
    walkers = sim_params.get('walkers', 100)
    trials = sim_params.get('trials', 10)  # Use this variable 'trials'
    dt = sim_params.get('dt', 0.0001)  # Use this variable 'dt'
    step_size = sim_params.get('step_size', 0.001)  # Use this variable 'step_size'

    # Grid Params
    grid_params = config.get('grid', {})
    grid_size = grid_params.get('size', 200)  # Use this variable 'grid_size'
    grid_min = grid_params.get('min', -1.0)
    grid_max = grid_params.get('max', 1.0)

    # Boundary Params
    boundary_params = config.get('boundaries', {})
    use_pbc = boundary_params.get('pbc', False)  # Use this variable 'use_pbc'
    check_bounds = boundary_params.get('check_bounds', False) if not use_pbc else False  # Use 'check_bounds'

    # CTRW and Alpha Params
    ctrw_config = config.get('ctrw', {})
    use_ctrw = ctrw_config.get('enabled', False)
    selected_alpha_func = None
    alpha_params = {}
    alpha_type = 'N/A'
    ctrw_config = config.get('ctrw', {})
    use_ctrw = ctrw_config.get('enabled', False)
    selected_alpha_func = None  # Will be populated below if use_ctrw is True
    alpha_params = {}
    alpha_type = 'N/A'  # Default value, will be updated below if use_ctrw is True

    # === CORRECTED BLOCK START ===
    if use_ctrw:  # Correctly check the boolean variable
        # Look for 'alpha_function' section as defined in your YAML
        alpha_function_section = ctrw_config.get('alpha_function', {})
        if not alpha_function_section:
            print("Warning: CTRW is enabled, but 'alpha_function' section is missing or empty in config.")
        else:
            # Get type and params from the 'alpha_function' section
            alpha_type = alpha_function_section.get('type', 'N/A')  # Get type, default to 'N/A' if missing
            alpha_params = alpha_function_section.get('params', {})

            # Validate the retrieved alpha_type
            if alpha_type not in AVAILABLE_ALPHA_FUNCTIONS:
                print(
                    f"Error: Unknown or missing alpha function type '{alpha_type}' in config. Available: {list(AVAILABLE_ALPHA_FUNCTIONS.keys())}")
                # Optionally set use_ctrw back to False or exit
                use_ctrw = False  # Safer to disable CTRW if type is invalid
                alpha_type = 'N/A'  # Reset type if invalid
                # import sys; sys.exit(1) # Or exit if preferred
            else:
                # If type is valid, store the function object (optional here, as it's passed later)
                selected_alpha_func = AVAILABLE_ALPHA_FUNCTIONS[alpha_type]
                print(f"  Successfully read alpha function type: {alpha_type}")  # Add confirmation
    # === CORRECTED BLOCK END ===


    # Disorder Params and Mode Determination
    disorder_config = config.get('disorder', {})
    disorder_type = disorder_config.get('type', 'none')
    disorder_params = disorder_config.get('params', {})
    selected_disorder_func = None
    disorder_mode = 'none'  # Default

    if disorder_type != 'none':
        if disorder_type in AVAILABLE_STANDARD_DISORDER_FUNCTIONS:
            selected_disorder_func = AVAILABLE_STANDARD_DISORDER_FUNCTIONS[disorder_type]
            disorder_mode = 'standard'
        elif disorder_type in AVAILABLE_DIRECTION_DISORDER_FUNCTIONS:
            selected_disorder_func = AVAILABLE_DIRECTION_DISORDER_FUNCTIONS[disorder_type]
            disorder_mode = 'direction_only'
            if use_ctrw:  # Disable CTRW if only directional disorder
                print("Warning: CTRW disabled because 'direction_only' disorder mode selected.")
                use_ctrw = False
        else:
            print(
                f"Error: Unknown disorder type '{disorder_type}'. Available: {list(ALL_AVAILABLE_DISORDER_FUNCTIONS.keys())}")
            import sys;

            sys.exit(1)


    # Animation Params
    anim_config = config.get('animation', {})
    run_animation = anim_config.get('enabled', False)  # Use this variable 'run_animation'
    save_animation = anim_config.get('save', False)
    anim_steps = anim_config.get('steps', 500)  # Use this variable 'anim_steps'
    anim_walker = anim_config.get('walker_index', 0)  # Use this variable 'anim_walker'
    anim_filename = anim_config.get('filename', 'walk_animation.gif')  # Use this variable 'anim_filename'

    # Histogram Params
    hist_config = config.get('histograms', {})
    run_histograms = hist_config.get('enabled', False)  # Use this variable 'run_histograms'
    hist_steps_to_plot = hist_config.get('steps_to_plot', [])  # Use 'hist_steps_to_plot'

    # --- Print Loaded Configuration Summary ---
    # (This part correctly uses the local variables)
    print("\n--- Simulation Configuration ---")
    print(f"  Steps: {steps}, Walkers: {walkers}, Trials: {trials}")
    print(f"  dt: {dt}, Step Size: {step_size}")
    print(f"  Grid: {grid_size}x{grid_size} from {grid_min} to {grid_max}")
    print(f"  Boundaries: PBC={use_pbc}, CheckBounds={check_bounds}")
    print(f"  Disorder: Mode='{disorder_mode}', Type='{disorder_type}', Params={disorder_params}")
    print(f"  CTRW: Enabled={use_ctrw}")
    if use_ctrw: print(f"    Alpha Function: Type='{alpha_type}', Params={alpha_params}")

    print("-" * 30)
    print(f"  Animation: Run={run_animation}, Save={save_animation}")
    print(f"  Histograms: Run={run_histograms}, Steps={hist_steps_to_plot}")
    print("-" * 30)

    # --- Setup Grid ---
    print(f"Setting up grid ({grid_size}x{grid_size})...")
    XV, YV = np.meshgrid(np.linspace(grid_min, grid_max, grid_size),
                         np.linspace(grid_min, grid_max, grid_size))
    print("Grid setup done.")

    # --- Run Parallel Simulation ---
    print(f"\n--- Running Parallel Simulation ({trials} Trials) ---")
    # *** Use local variables loaded from config, NOT args.***
    avg_msd, time_axis = main_parallel(
        num_trials_total=trials, num_steps=steps, num_walkers=walkers, step=step_size, dt=dt,
        disorder_type=disorder_type,  # Pass type name
        disorder_params=disorder_params,
        alpha_type=alpha_type,  # Pass type name
        alpha_params=alpha_params,
        disorder_mode=disorder_mode,
        xv=XV, yv=YV, use_ctrw_flag=use_ctrw,
        use_pbc_flag=use_pbc, check_bounds_flag=check_bounds
    )

    # --- Quantitative Analysis ---
    if avg_msd is not None and time_axis is not None:
        print("\n" + "=" * 30)
        print(" Quantitative Analysis Results")
        print("=" * 30)

        # Define fit range (e.g., last half of the data, avoiding first few points)
        # Ensure steps is defined correctly from your config loading
        min_fit_step = max(10, int(steps // 3))  # Use int() for safety if steps is float
        max_fit_step = int(steps)
        print(f"Analysis Range Steps: [{min_fit_step}, {max_fit_step}]")  # Print range once

        # Initialize fit results to NaN
        fitted_alpha_msd = np.nan
        fit_intercept_msd = np.nan
        fitted_alpha_msd_t = np.nan

        # Check if the range is valid
        if max_fit_step > min_fit_step and len(time_axis) > max_fit_step:
            # Get indices corresponding to the step range
            idx_min = min_fit_step
            idx_max = max_fit_step  # Use index corresponding to max_fit_step
            time_fit_range = time_axis[idx_min: idx_max + 1]
            msd_fit_range = avg_msd[idx_min: idx_max + 1]

            # --- 1. Fit Log-Log MSD to find alpha exponent (Overall Fit) ---
            print("\n--- Method 1: Log-Log MSD Fit (log(MSD) vs log(t)) ---")
            valid_fit_indices_msd = (time_fit_range > 1e-15) & (msd_fit_range > 1e-15)
            if np.sum(valid_fit_indices_msd) >= 2:
                log_time_msd = np.log(time_fit_range[valid_fit_indices_msd])
                log_msd = np.log(msd_fit_range[valid_fit_indices_msd])
                try:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(log_time_msd, log_msd)
                    fitted_alpha_msd = slope  # Store overall alpha
                    fit_intercept_msd = intercept  # Store overall intercept (log(C))
                    print(f"  Estimated Overall Alpha (Slope) = {fitted_alpha_msd:.4f}")
                    print(f"  Standard Error                  = {std_err:.4f}")
                    print(f"  R-squared                       = {r_value ** 2:.4f}")
                except Exception as e:
                    print(f"  Error during log-log MSD fit: {e}")
                    fitted_alpha_msd = np.nan  # Ensure NaN on error
                    fit_intercept_msd = np.nan
            else:
                print(f"  Not enough valid data points ({np.sum(valid_fit_indices_msd)}) for log-log MSD fit.")

            # --- 2. Fit Log-Log MSD/Time to find alpha-1 exponent ---
            print("\n--- Method 2: Log-Log MSD/Time Fit (log(MSD/t) vs log(t)) ---")
            # Avoid division by zero in MSD/t calculation
            valid_time_for_div = (time_fit_range > 1e-15)
            if np.any(valid_time_for_div):
                msd_over_time_fit_range = msd_fit_range[valid_time_for_div] / time_fit_range[valid_time_for_div]
                time_for_msd_t_fit = time_fit_range[valid_time_for_div]
                # Also check MSD/t > 0 for log
                valid_fit_indices_msd_t = (msd_over_time_fit_range > 1e-15)

                if np.sum(valid_fit_indices_msd_t) >= 2:
                    log_time_msd_t = np.log(time_for_msd_t_fit[valid_fit_indices_msd_t])
                    log_msd_over_time = np.log(msd_over_time_fit_range[valid_fit_indices_msd_t])
                    try:
                        slope_alpha_minus_1, intercept_b, r_value_b, p_value_b, std_err_b = stats.linregress(
                            log_time_msd_t,
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
                print(f"  Not enough valid time points > 0 in range for MSD/t calculation.")

        else:
            print(
                f"  Fit range steps [{min_fit_step}, {max_fit_step}] invalid or insufficient data length ({len(time_axis)} points).")

        # --- 3. Calculate Effective Diffusion Coefficient ---
        print("\n--- Effective Diffusion Coefficient (D_eff = MSD / 4t) ---")
        # Recalculate D_eff over the same fit range used for alpha
        if max_fit_step > min_fit_step and len(time_axis) > max_fit_step:
            # Use idx_min, idx_max defined earlier
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

                print(f"  Analysis Range Steps: [{min_fit_step}, {max_fit_step}]")
                print(f"  Average D_eff in range = {avg_d_eff:.4e}")
                print(f"  Final D_eff in range   = {final_d_eff:.4e}")
            else:
                print("  No valid time points > 0 in range for D_eff calculation.")
                avg_d_eff = np.nan
                final_d_eff = np.nan
        else:
            # This message might be redundant if the outer check already printed
            # print(f"  Fit range [{min_fit_step}, {max_fit_step}] invalid or insufficient data for D_eff.")
            avg_d_eff = np.nan
            final_d_eff = np.nan

        # --- 4. Theoretical Comparison (for ordered case) ---
        # Ensure step_size and dt are defined from config loading
        D_theory_ordered = step_size ** 2 / (4 * dt)
        print("\n--- Theoretical Comparison ---")
        print(f"  Theoretical D (Ordered Walk) = {D_theory_ordered:.4e}")
        # Ensure disorder_type and use_ctrw are defined from config loading
        if disorder_type == 'none' and not use_ctrw:
            print(f"  (Simulation matches theoretical D if Avg D_eff -> Theoretical D)")
        else:
            print(f"  (Disorder/CTRW expected to reduce D_eff compared to theoretical)")

        print("=" * 30 + "\n")

        # --- Plotting ---
        print("--- Plotting Averaged Results ---")

        # Plot MSD/Time
        plt.figure(figsize=(10, 6))
        valid_div_plot = time_axis > 1e-15
        if np.any(valid_div_plot):
            msd_over_time_plot = avg_msd[valid_div_plot] / time_axis[valid_div_plot]
            # Ensure alpha_type is defined from config loading
            plt.plot(time_axis[valid_div_plot], msd_over_time_plot,
                     label=f'MSD/Time ({disorder_type}, alpha={alpha_type})')
            plt.axhline(4 * D_theory_ordered, color='r', linestyle='--', alpha=0.7,
                        label=f'4 * D_theory (Ordered) = {4 * D_theory_ordered:.2e}')
        plt.xlabel('Time (s)')
        plt.ylabel('MSD / Time')
        plt.title(f'Avg Effective Diffusion Coefficient ({trials} Trials)')  # Ensure trials is defined
        plt.grid(True)
        plt.legend()
        plt.show()

        # Plot Log-Log MSD
        plt.figure(figsize=(10, 6))
        # Filter for valid log values right at the start
        valid_log_plot = (time_axis > 1e-15) & (avg_msd > 1e-15)
        time_axis_valid_plot = time_axis[valid_log_plot]
        avg_msd_valid_plot = avg_msd[valid_log_plot]
        log_time_full_plot = np.log(time_axis_valid_plot)
        log_msd_full_plot = np.log(avg_msd_valid_plot)

        if np.any(valid_log_plot):
            n_points_plot = len(time_axis_valid_plot)  # Number of valid points for plotting

            # Plot the actual MSD data first
            plt.loglog(time_axis_valid_plot, avg_msd_valid_plot, 'o', markersize=3, alpha=0.6,
                       label=f'MSD Data ({disorder_type}, alpha={alpha_type})')

            # Plot theoretical line (Slope=1)
            slope_1_line = 4 * D_theory_ordered * time_axis_valid_plot
            plt.loglog(time_axis_valid_plot, slope_1_line, 'r--', alpha=0.7,
                       label=f'Slope=1 (Theory D={D_theory_ordered:.2e})')

            # --- Conditional Fitting for Plotting---
            # Ensure use_ctrw and alpha_type are defined from config loading
            if use_ctrw and alpha_type in ['gaussian', 'linear_gradient']:
                print(f"--- Plotting Two Log-Log Fits for {alpha_type} alpha ---")

                # --- Define Plot Fit Ranges (Based on valid plot points) ---
                # Short time: e.g., from index 5 to 10% of points (min 10 points)
                start_short_plot = 0
                end_short_plot = max(start_short_plot + 9, n_points_plot // 5000)  # Ensure at least 10 points if possible
                # Long time: e.g., from 50% of points to the end (min 10 points)
                start_long_plot = max(n_points_plot // 100, end_short_plot + 1)  # Ensure no overlap
                end_long_plot = n_points_plot - 1

                # --- Fit 1: Short Time ---
                # Check if range is valid within the plotted data
                if start_short_plot <= end_short_plot and (
                        end_short_plot - start_short_plot) >= 1:  # Need at least 2 points for fit
                    try:
                        slope_short, intercept_short, r_short, p_short, stderr_short = stats.linregress(
                            log_time_full_plot[start_short_plot:end_short_plot + 1],
                            log_msd_full_plot[start_short_plot:end_short_plot + 1]
                        )
                        fitted_alpha_short = slope_short
                        # Use the valid time axis subset for plotting the line
                        time_short_range_plot = time_axis_valid_plot[start_short_plot:end_short_plot + 1]
                        fit_line_short = np.exp(intercept_short) * (time_short_range_plot ** fitted_alpha_short)
                        plt.loglog(time_short_range_plot, fit_line_short, 'g-', linewidth=2,
                                   label=f'Short Time Fit (α={fitted_alpha_short:.3f})')
                        print(f"  Plotting Short Time Fit (indices {start_short_plot}-{end_short_plot})")
                    except Exception as e:
                        print(f"  Error during short time fit for plot: {e}")
                else:
                    print(
                        f"  Not enough points for short time fit plot (Range indices: {start_short_plot}-{end_short_plot}, Available: {n_points_plot})")

                # --- Fit 2: Long Time ---
                # Check if range is valid within the plotted data
                if start_long_plot <= end_long_plot and (
                        end_long_plot - start_long_plot) >= 1:  # Need at least 2 points for fit
                    try:
                        slope_long, intercept_long, r_long, p_long, stderr_long = stats.linregress(
                            log_time_full_plot[start_long_plot:end_long_plot + 1],
                            log_msd_full_plot[start_long_plot:end_long_plot + 1]
                        )
                        fitted_alpha_long = slope_long
                        # Use the valid time axis subset for plotting the line
                        time_long_range_plot = time_axis_valid_plot[start_long_plot:end_long_plot + 1]
                        fit_line_long = np.exp(intercept_long) * (time_long_range_plot ** fitted_alpha_long)
                        plt.loglog(time_long_range_plot, fit_line_long, 'm--', linewidth=2,
                                   label=f'Long Time Fit (α={fitted_alpha_long:.3f})')
                        print(f"  Plotting Long Time Fit (indices {start_long_plot}-{end_long_plot})")
                    except Exception as e:
                        print(f"  Error during long time fit for plot: {e}")
                else:
                    print(
                        f"  Not enough points for long time fit plot (Range indices: {start_long_plot}-{end_long_plot}, Available: {n_points_plot})")

            else:
                # --- Original Single Fit Plotting Logic ---
                print("--- Plotting Single Overall Log-Log Fit ---")
                # Use the overall fit results calculated in the analysis section
                if not np.isnan(fitted_alpha_msd) and not np.isnan(fit_intercept_msd):
                    # Calculate fitted line over the whole valid range
                    fit_line_msd_plot = np.exp(fit_intercept_msd) * (time_axis_valid_plot ** fitted_alpha_msd)
                    plt.loglog(time_axis_valid_plot, fit_line_msd_plot, 'g:', alpha=0.9, linewidth=2,
                               label=f'Overall Fit (α={fitted_alpha_msd:.3f})')
                    print(f"  Plotting Overall Fit (alpha = {fitted_alpha_msd:.4f})")
                else:
                    print("  Overall fit calculation failed or skipped, not plotted.")

            # --- Final Plot Settings ---
            plt.xlabel('Time (s)')
            plt.ylabel('MSD')
            plt.title(f'Avg Mean Squared Displacement (Log-Log, {trials} Trials)')  # Ensure trials defined
            plt.grid(True, which='both');
            plt.legend()
            plt.show()

        else:
            print("No valid data points for log-log plotting.")

    else:
        print("Parallel simulation failed or produced no results, skipping analysis and plotting.")

    # --- End of the main analysis and plotting block ---
    # --- Optional: Run Single Trial for Animation ---
    # *** Use local variables loaded from config ***
    if run_animation:
        print("\n--- Running Single Trial for Animation ---")
        rw_anim = RandomWalk(
            num_steps=anim_steps,  # Use 'anim_steps'
            num_walkers=walkers,
            step=step_size, dt=dt, xv=XV, yv=YV,
            disorder_function=selected_disorder_func,
            disorder_params=disorder_params,alpha_function=selected_alpha_func,alpha_params=alpha_params,disorder_mode=disorder_mode,
            use_ctrw=use_ctrw,
            use_pbc=use_pbc,
            check_bounds=check_bounds,
            store_history=True
        )
        use_disorder_anim = (selected_disorder_func is not None)
        rw_anim.trajectories()
        rw_anim.animate_trajectory(
            walker_index=anim_walker,  # Use 'anim_walker'
            save_animation=save_animation,  # Use 'save_animation'
            filename=anim_filename  # Use 'anim_filename'
        )

    # --- Optional: Run Single Trial for Histograms ---
    # *** Use local variables loaded from config ***
    if run_histograms and hist_steps_to_plot:
        print("\n--- Running Single Trial for Histograms ---")
        max_hist_step = max(hist_steps_to_plot) if hist_steps_to_plot else 0
        # Ensure simulation runs long enough for both main analysis AND histograms
        hist_run_steps = max(steps, max_hist_step)

        # Ensure 'alpha_type' is defined in this scope, loaded from your config.
        # It was defined in your main section:
        # if use_ctrw:
        #     alpha_function_section = ctrw_config.get('alpha_function', {})
        #     alpha_type = alpha_function_section.get('type', 'N/A')
        #     alpha_params = alpha_function_section.get('params', {})
        # else:
        #     alpha_type = 'N/A' # Or some default if use_ctrw is false
        #     alpha_params = {}

        rw_hist = RandomWalk(
            num_steps=hist_run_steps,
            num_walkers=walkers,
            step=step_size,
            dt=dt,
            xv=XV,
            yv=YV,
            disorder_function=selected_disorder_func,
            disorder_params=disorder_params,
            # REMOVE: alpha_function=selected_alpha_func,
            alpha_type=alpha_type,  # CORRECT: Pass the string 'alpha_type'
            alpha_params=alpha_params,
            disorder_mode=disorder_mode,
            use_ctrw=use_ctrw,
            use_pbc=use_pbc,
            check_bounds=check_bounds,
            store_history=True
        )
        # use_disorder_hist = (selected_disorder_func is not None) # This line seems okay
        rw_hist.trajectories()
        print("Generating Histograms...")
        for step_to_plot in hist_steps_to_plot:  # Use 'hist_steps_to_plot'
            if step_to_plot <= hist_run_steps:
                rw_hist.plot_position_histograms(time_step=step_to_plot)
            else:
                print(f"Warning: Requested hist step {step_to_plot} > sim length {hist_run_steps}")
        plt.show()

    print("\nSimulation Finished.")