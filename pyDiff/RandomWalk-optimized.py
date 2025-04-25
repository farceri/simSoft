import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numba
import multiprocessing # Import the module
import time as timer # To time the execution
import os # To potentially get CPU count

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




IDX_P_X, IDX_M_X, IDX_P_Y, IDX_M_Y, IDX_REST = 0, 1, 2, 3, 4
NUM_DIRECTIONS = 4 # Number of movement directions

def my_spatial_disorder_vectorized(x, y, rest_fixed=0.9, **kwargs):
    """
    Vectorized version of my_spatial_disorder (using the active part).
    Applies a fixed uniform resting probability across the grid.

    Args:
        x (np.ndarray): Meshgrid of x-coordinates (shape ny, nx).
        y (np.ndarray): Meshgrid of y-coordinates (shape ny, nx).
        rest_fixed (float): The fixed resting probability to apply everywhere.
        **kwargs: Catches unused parameters passed during precomputation.

    Returns:
        np.ndarray: Array of probabilities (shape ny, nx, 5) with dtype float32.
    """
    # Ensure rest_fixed is valid
    rest_prob_val = np.clip(rest_fixed, 0.0, 1.0)

    # Create an array of the same shape as x (or y) filled with the rest probability
    rest_prob = np.full(x.shape, rest_prob_val, dtype=np.float32)

    # Calculate movement probability (element-wise)
    move_prob_total = 1.0 - rest_prob
    # Use np.maximum for vectorized max(0.0, ...)
    move_prob_each = np.maximum(0.0, move_prob_total / NUM_DIRECTIONS)

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs

def gaussian_rest_prob_vectorized(x, y, sigma=0.1, max_rest_strength=1.0, **kwargs):
    """
    Vectorized Gaussian resting probability centered at the origin.
    Ensures the returned probabilities always sum to 1.0.

    Args:
        x (np.ndarray): Meshgrid of x-coordinates (shape ny, nx).
        y (np.ndarray): Meshgrid of y-coordinates (shape ny, nx).
        sigma (float): Standard deviation of the Gaussian distribution.
        max_rest_strength (float): The maximum resting probability at the origin (must be <= 1).
        **kwargs: Catches unused parameters passed during precomputation.


    Returns:
        np.ndarray: Array of probabilities (shape ny, nx, 5) with dtype float32.
    """
    # Ensure max_rest_strength is valid
    max_rest_strength = min(max_rest_strength, 1.0) # Cannot be more than 1

    # Calculate the resting probability based on Gaussian decay (element-wise)
    raw_rest_prob = max_rest_strength * np.exp(-(x**2 + y**2) / (2 * sigma**2))

    # Ensure rest_prob is strictly within [0, 1]
    rest_prob = np.clip(raw_rest_prob, 0.0, 1.0).astype(np.float32)

    # Calculate movement probability (element-wise)
    move_prob_total = 1.0 - rest_prob
    move_prob_each = np.maximum(0.0, move_prob_total / NUM_DIRECTIONS)

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs


def uniform_rest_prob_vectorized(x, y, rest_level=0.9, **kwargs):
    """
    Vectorized uniform resting probability across the grid.
    Ensures the returned probabilities always sum to 1.0.

    Args:
        x (np.ndarray): Meshgrid of x-coordinates (shape ny, nx).
        y (np.ndarray): Meshgrid of y-coordinates (shape ny, nx).
        rest_level (float): The uniform resting probability level (clipped to [0, 1]).
        **kwargs: Catches unused parameters passed during precomputation.

    Returns:
        np.ndarray: Array of probabilities (shape ny, nx, 5) with dtype float32.
    """
    # Ensure rest_level is valid
    rest_prob_val = np.clip(rest_level, 0.0, 1.0)

    # Create an array filled with the rest probability
    rest_prob = np.full(x.shape, rest_prob_val, dtype=np.float32)

    # Calculate movement probability
    move_prob_total = 1.0 - rest_prob
    move_prob_each = np.maximum(0.0, move_prob_total / NUM_DIRECTIONS)

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs

def plateau_rest_prob_vectorized(x, y, plateau_radius=0.001, max_rest=0.5, decay_rate=5.0, **kwargs):
    """
    Vectorized version: Resting probability has a plateau near the origin.
    Ensures the returned probabilities always sum to 1.0.

     Args:
        x (np.ndarray): Meshgrid of x-coordinates (shape ny, nx).
        y (np.ndarray): Meshgrid of y-coordinates (shape ny, nx).
        plateau_radius (float): Radius of the central plateau.
        max_rest (float): Resting probability within the plateau.
        decay_rate (float): Exponential decay rate outside the plateau.
        **kwargs: Catches unused parameters passed during precomputation.

    Returns:
        np.ndarray: Array of probabilities (shape ny, nx, 5) with dtype float32.
    """
    distance_from_origin = np.sqrt(x**2 + y**2)

    # Calculate raw rest probability using np.where for conditional logic on arrays
    raw_rest_prob = np.where(
        distance_from_origin <= plateau_radius,
        max_rest, # Value if condition is true
        max_rest * np.exp(-decay_rate * (distance_from_origin - plateau_radius)) # Value if false
    )

    # Ensure rest_prob is in [0, 1]
    rest_prob = np.clip(raw_rest_prob, 0.0, 1.0).astype(np.float32)

    # Calculate movement probability (element-wise)
    move_prob_total = 1.0 - rest_prob
    move_prob_each = np.maximum(0.0, move_prob_total / NUM_DIRECTIONS)

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs

def multi_center_rest_prob_vectorized(x, y, center1=(0.001, 0.001), center2=(-0.001, -0.001),
                                      strength1=0.5, strength2=0.5, decay_rate=5.0, max_total_rest=0.9, **kwargs):
    """
    Vectorized version: High resting probability around multiple centers.
    Ensures the returned probabilities always sum to 1.0.

    Args:
        x (np.ndarray): Meshgrid of x-coordinates (shape ny, nx).
        y (np.ndarray): Meshgrid of y-coordinates (shape ny, nx).
        center1 (tuple): Coordinates (x, y) of the first center.
        center2 (tuple): Coordinates (x, y) of the second center.
        strength1 (float): Max strength of the first center's rest probability.
        strength2 (float): Max strength of the second center's rest probability.
        decay_rate (float): Gaussian decay rate for both centers.
        max_total_rest (float): Maximum allowed combined resting probability.
       **kwargs: Catches unused parameters passed during precomputation.

    Returns:
        np.ndarray: Array of probabilities (shape ny, nx, 5) with dtype float32.
    """
    # Calculate distance squared to each center
    dist_sq1 = (x - center1[0])**2 + (y - center1[1])**2
    dist_sq2 = (x - center2[0])**2 + (y - center2[1])**2

    # Calculate rest probability contribution from each center
    center1_rest = strength1 * np.exp(-decay_rate * dist_sq1)
    center2_rest = strength2 * np.exp(-decay_rate * dist_sq2)

    # Combine contributions and clip to the overall maximum allowed rest probability
    raw_rest_prob = np.clip(center1_rest + center2_rest, 0.0, max_total_rest)

    # Ensure rest_prob is in [0, 1] (redundant if max_total_rest <= 1, but safe)
    rest_prob = np.clip(raw_rest_prob, 0.0, 1.0).astype(np.float32)

    # Calculate movement probability (element-wise)
    move_prob_total = 1.0 - rest_prob
    move_prob_each = np.maximum(0.0, move_prob_total / NUM_DIRECTIONS)

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs


def exponential_rest_prob_vectorized(x, y, decay_rate=0.1, max_rest=0.95, **kwargs):
    """
    Vectorized version: Resting probability decreases exponentially from the origin.
    Ensures the returned probabilities always sum to 1.0 using EQUAL move probabilities.

    Args:
        x (np.ndarray): Meshgrid of x-coordinates (shape ny, nx).
        y (np.ndarray): Meshgrid of y-coordinates (shape ny, nx).
        decay_rate (float): Exponential decay rate based on distance.
        max_rest (float): Maximum resting probability at the origin.
        **kwargs: Catches unused parameters passed during precomputation.

    Returns:
        np.ndarray: Array of probabilities (shape ny, nx, 5) with dtype float32.
    """
    distance_from_origin = np.sqrt(x**2 + y**2)

    # Calculate the raw resting probability
    raw_rest_prob = max_rest * np.exp(-decay_rate * distance_from_origin)

    # Ensure rest_prob is in [0, 1]
    rest_prob = np.clip(raw_rest_prob, 0.0, 1.0).astype(np.float32)

    # Calculate movement probability (element-wise) - distributes remaining probability equally
    move_prob_total = 1.0 - rest_prob
    move_prob_each = np.maximum(0.0, move_prob_total / NUM_DIRECTIONS)

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs

def boundary_dependent_rest_prob_vectorized(x, y, x_bounds=(-1.0, 1.0), y_bounds=(-1.0, 1.0),
                                            boundary_strength=0.5, decay_rate=5.0, max_total_rest=0.9, **kwargs):
    """
    Vectorized version: Resting probability increases near the boundaries.
    Ensures the returned probabilities always sum to 1.0.

    Args:
        x (np.ndarray): Meshgrid of x-coordinates (shape ny, nx).
        y (np.ndarray): Meshgrid of y-coordinates (shape ny, nx).
        x_bounds (tuple): (min_x, max_x) defining the boundary region.
        y_bounds (tuple): (min_y, max_y) defining the boundary region.
        boundary_strength (float): Strength scaling factor for boundary effect.
        decay_rate (float): Exponential decay rate from the boundary.
        max_total_rest (float): Maximum allowed combined resting probability.
        **kwargs: Catches unused parameters passed during precomputation.

    Returns:
        np.ndarray: Array of probabilities (shape ny, nx, 5) with dtype float32.
    """
    # Calculate distance from boundaries
    dist_to_min_x = np.abs(x - x_bounds[0])
    dist_to_max_x = np.abs(x - x_bounds[1])
    dist_to_min_y = np.abs(y - y_bounds[0])
    dist_to_max_y = np.abs(y - y_bounds[1])

    # Calculate rest probability contribution from each boundary edge
    rest_prob_min_x = boundary_strength * np.exp(-decay_rate * dist_to_min_x)
    rest_prob_max_x = boundary_strength * np.exp(-decay_rate * dist_to_max_x)
    rest_prob_min_y = boundary_strength * np.exp(-decay_rate * dist_to_min_y)
    rest_prob_max_y = boundary_strength * np.exp(-decay_rate * dist_to_max_y)

    # Combine contributions (summing effect from all boundaries)
    # Clip to the overall maximum allowed rest probability
    raw_rest_prob = np.clip(rest_prob_min_x + rest_prob_max_x + rest_prob_min_y + rest_prob_max_y, 0.0, max_total_rest)

    # Ensure rest_prob is in [0, 1]
    rest_prob = np.clip(raw_rest_prob, 0.0, 1.0).astype(np.float32)

    # Calculate movement probability (element-wise)
    move_prob_total = 1.0 - rest_prob
    move_prob_each = np.maximum(0.0, move_prob_total / NUM_DIRECTIONS)

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs









def corrected_gaussian_rest_prob_vectorized(x, y, sigma=1.0, max_rest_strength=0.95, **kwargs):
    """
    Vectorized version: Calculates probabilities for all input x, y points.
    x, y are expected to be NumPy arrays (like meshgrids xv, yv).
    Returns an array of shape (*x.shape, 5).
    CORRECTED: Uses np.minimum instead of Python min.
    """
    # Use np.minimum for compatibility with potential array operations
    max_rest_strength = np.minimum(max_rest_strength, 1.0)

    # These operations work element-wise on arrays x, y
    rest_prob = max_rest_strength * np.exp(-(x**2 + y**2) / (2 * sigma**2))
    rest_prob = np.clip(rest_prob, 0.0, 1.0)

    move_prob_total = 1.0 - rest_prob
    move_prob_each = np.maximum(0.0, move_prob_total / 4.0) # np.maximum is correct

    # Create the output array (ny, nx, 5)
    out_shape = x.shape + (5,)
    probs = np.zeros(out_shape, dtype=np.float32)

    probs[..., IDX_P_X] = move_prob_each
    probs[..., IDX_M_X] = move_prob_each
    probs[..., IDX_P_Y] = move_prob_each
    probs[..., IDX_M_Y] = move_prob_each
    probs[..., IDX_REST] = rest_prob

    return probs




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
    time_axis = np.arange(num_steps + 1)
    return avg_msd, time_axis
# =============================================================================






"""
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
        rw_gaussian_rest_new = RandomWalk(disorder_function=corrected_gaussian_rest_prob_vectorized, num_steps=num_steps)
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
"""

if __name__ == "__main__":

    # --- Simulation Parameters ---
    ENABLE_CTRW = True  # <<< SET TO True TO ENABLE CTRW, False FOR STANDARD REST >>>
    CTRW_ALPHA = 0.7  # <<< SET desired exponent if ENABLE_CTRW is True >>>




    # --- Simulation Parameters ---
    NUM_TRIALS = 20 # Number of parallel trials
    NUM_STEPS = 10000 # Number of steps per trial
    NUM_WALKERS = 100
    STEP_SIZE = 0.001
    TIME_STEP_DT = 0.0001

    # Define grid (consider making it smaller if memory/precomputation is slow)
    grid_size = 2000 # Example: Use a smaller grid for testing
    print(f"Setting up grid ({grid_size}x{grid_size})...")
    XV, YV = np.meshgrid(np.linspace(-1, 1, grid_size), np.linspace(-1, 1, grid_size))
    print("Grid setup done.")

    # --- Select Disorder Function and Parameters ---
    # Example: Corrected Gaussian
    DISORDER_FUNC = corrected_gaussian_rest_prob_vectorized # Use the vectorized version
    DISORDER_PARAMS = {'sigma': 1.0, 'max_rest_strength': 0.95}
    if ENABLE_CTRW:
        DISORDER_PARAMS['ctrw_alpha'] = CTRW_ALPHA # Add alpha only if CTRW is on



    # Example: No disorder
    #DISORDER_FUNC = None
    #DISORDER_PARAMS = {}
    #if ENABLE_CTRW:
     #   DISORDER_PARAMS['ctrw_alpha'] = CTRW_ALPHA # Add alpha only if CTRW is on


    # --- Run the parallel simulation ---
    avg_msd, time_axis = main_parallel(
        num_trials_total=NUM_TRIALS,
        num_steps=NUM_STEPS,
        num_walkers=NUM_WALKERS,
        step=STEP_SIZE,
        dt=TIME_STEP_DT,
        disorder_function=DISORDER_FUNC,
        disorder_params=DISORDER_PARAMS,
        xv=XV,
        yv=YV,
        use_ctrw_flag=ENABLE_CTRW  # Pass the flag
    )

    # --- Quantitative Analysis & Plotting Results ---
    if avg_msd is not None and time_axis is not None:
        print("\n" + "=" * 30)
        print(" Quantitative Analysis Results")
        print("=" * 30)

        # === Method 1: Fit Log-Log MSD ===
        print("\n--- Analysis Method 1: Log-Log MSD Fit ---")
        N_min_fit = NUM_STEPS // 2  # Example: Fit last half
        N_max_fit = NUM_STEPS
        N_min_fit = max(100, N_min_fit)  # Ensure minimum range, avoid early transients

        if N_max_fit <= N_min_fit:
            print("Not enough time steps for fitting range.")
        else:
            time_fit_range = time_axis[N_min_fit: N_max_fit + 1]
            msd_fit_range = avg_msd[N_min_fit: N_max_fit + 1]
            valid_indices = (time_fit_range > 0) & (msd_fit_range > 0)
            if np.sum(valid_indices) < 2:
                print("Not enough valid data points for log-log fit.")
            else:
                time_log = np.log(time_fit_range[valid_indices])
                msd_log = np.log(msd_fit_range[valid_indices])
                slope, intercept, r_value, p_value, std_err = stats.linregress(time_log, msd_log)
                alpha_estimate = slope
                print(f"Fit Range N = [{N_min_fit}, {N_max_fit}]")
                print(f"  Estimated alpha (slope) = {alpha_estimate:.4f}")
                print(f"  Standard Error          = {std_err:.4f}")
                print(f"  R-squared               = {r_value ** 2:.4f}")

        # === Method 2: Fit Log-Log MSD/N ===
        print("\n--- Analysis Method 2: Log-Log MSD/N Fit ---")
        time_eff = time_axis[1:]
        msd_over_n = avg_msd[1:] / time_eff

        # Use same fit range N_min_fit, N_max_fit
        if N_max_fit <= N_min_fit:
            print("Not enough time steps for fitting range.")
        else:
            idx_min = N_min_fit - 1;
            idx_max = N_max_fit - 1
            time_fit_range_eff = time_eff[idx_min: idx_max + 1]
            msd_over_n_fit_range = msd_over_n[idx_min: idx_max + 1]
            valid_indices_eff = (time_fit_range_eff > 0) & (msd_over_n_fit_range > 0)
            if np.sum(valid_indices_eff) < 2:
                print("Not enough valid data points for log-log MSD/N fit.")
            else:
                time_log_eff = np.log(time_fit_range_eff[valid_indices_eff])
                msd_over_n_log = np.log(msd_over_n_fit_range[valid_indices_eff])
                slope_b, intercept_c, r_value_b, p_value_b, std_err_b = stats.linregress(time_log_eff, msd_over_n_log)
                alpha_minus_1_estimate = slope_b
                alpha_estimate_from_b = slope_b + 1.0
                print(f"Fit Range N = [{N_min_fit}, {N_max_fit}]")
                print(f"  Estimated alpha-1 (slope) = {alpha_minus_1_estimate:.4f}")
                print(f"  Implied alpha             = {alpha_estimate_from_b:.4f}")
                print(f"  Standard Error (of slope) = {std_err_b:.4f}")
                print(f"  R-squared                 = {r_value_b ** 2:.4f}")
        print("=" * 30 + "\n")

        # --- Your existing plotting code ---
        print("Plotting results...")
        mode_label = f"CTRW alpha={CTRW_ALPHA}" if ENABLE_CTRW else "Standard Rest"
        func_name = DISORDER_FUNC.__name__ if DISORDER_FUNC else "No Resting"

        # Plot MSD/Time
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(time_axis[1:], avg_msd[1:] / time_axis[1:], label=f'{func_name} ({mode_label})')
        ax1.set_xlabel('Time Step')
        ax1.set_ylabel('MSD / Time Step')
        ax1.set_title(f'Avg Effective Diffusion Coefficient ({NUM_TRIALS} Trials)')
        ax1.grid(True);
        ax1.legend()
        plt.savefig(f"msd_over_time_{'ctrw' if ENABLE_CTRW else 'std'}.png")

        # Plot Log-Log MSD
        fig2, ax2 = plt.subplots(figsize=(10, 6))
        valid = (time_axis > 0) & (avg_msd > 0)
        ax2.loglog(time_axis[valid], avg_msd[valid], label=f'{func_name} ({mode_label})')
        # Add slope=1 line for reference
        if np.any(valid):
            first_msd = avg_msd[valid][0];
            first_time = time_axis[valid][0]
            # Use calculated alpha if available and reliable, otherwise default to 1 for guide
            guide_alpha = alpha_estimate if 'alpha_estimate' in locals() and not np.isnan(alpha_estimate) else 1.0
            # Generate fit line using results from Method 1
            if 'alpha_estimate' in locals() and not np.isnan(alpha_estimate):
                fit_line = np.exp(intercept) * (time_axis[valid] ** alpha_estimate)
                ax2.loglog(time_axis[valid], fit_line, 'r--', alpha=0.7, label=f'Fit (alpha={alpha_estimate:.3f})')
            else:  # Fallback slope=1 guide if fit failed
                slope_1_line = (first_msd / first_time) * time_axis[valid]
                ax2.loglog(time_axis[valid], slope_1_line, 'r--', alpha=0.7, label='Slope=1 guide')

            slope_1_line = (first_msd / first_time) * time_axis[valid]
            ax2.loglog(time_axis[valid], slope_1_line, 'g--', alpha=0.7, label='Slope=1 standard diffusion')
        ax2.set_xlabel('Time Step');
        ax2.set_ylabel('MSD')
        ax2.set_title(f'Avg Mean Squared Displacement (Log-Log, {NUM_TRIALS} Trials)')
        ax2.grid(True, which='both');
        ax2.legend()
        plt.savefig(f"msd_loglog_{'ctrw' if ENABLE_CTRW else 'std'}.png")

        plt.show()
    else:
        print("Simulation failed, skipping analysis and plotting.")
    """

    print("\nRunning a single trial for animation...")


    # Use the same simulation parameters
    ENABLE_CTRW_ANIM = True  # Or True, depending on what you want to animate
    CTRW_ALPHA_ANIM = 0.7
    NUM_STEPS_ANIM = 1000  # Can be shorter for faster animation generation
    NUM_WALKERS_ANIM = 100
    STEP_SIZE_ANIM = 0.001
    TIME_STEP_DT_ANIM = 0.0001
    grid_size_anim = 200  # Smaller grid might be okay for visualization
    XV_anim, YV_anim = np.meshgrid(np.linspace(-1, 1, grid_size_anim), np.linspace(-1, 1, grid_size_anim))

    # Select disorder function for animation run
    #DISORDER_FUNC_ANIM = None  # Example: Ordered walk
    #DISORDER_PARAMS_ANIM = {}
    # OR
    DISORDER_FUNC_ANIM = corrected_gaussian_rest_prob_vectorized
    DISORDER_PARAMS_ANIM = {'sigma': 1.0, 'max_rest_strength': 0.95}
    if ENABLE_CTRW_ANIM:
        DISORDER_PARAMS_ANIM['ctrw_alpha'] = CTRW_ALPHA_ANIM

    # **** Crucial Step: Instantiate with store_history=True ****
    rw_anim = RandomWalk(
        num_steps=NUM_STEPS_ANIM,
        num_walkers=NUM_WALKERS_ANIM,
        step=STEP_SIZE_ANIM,
        dt=TIME_STEP_DT_ANIM,
        xv=XV_anim,
        yv=YV_anim,
        disorder_function=DISORDER_FUNC_ANIM,
        disorder_params=DISORDER_PARAMS_ANIM,
        use_ctrw=ENABLE_CTRW_ANIM,
        store_history=True  # <--- SET THIS TO TRUE
    )

    # Determine if disorder is used for this specific run
    use_disorder_anim = (DISORDER_FUNC_ANIM is not None)

    # Run the trajectories method for this single instance
    rw_anim.trajectories(use_disorder=use_disorder_anim)

    # --- Now you can animate ---
    print("Generating animation...")
    # Improve the check in animate_trajectory itself:
    if hasattr(rw_anim, 'all_positions') and rw_anim.all_positions is not None and len(rw_anim.all_positions) > 0:
        rw_anim.animate_trajectory(walker_index=0, interval=50, save_animation=False, filename="walk_animation.gif")
        print("Animation finished")
    else:
        print("Could not generate animation: Position history was not stored.")
    """