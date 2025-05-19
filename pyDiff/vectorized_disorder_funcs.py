import numpy as np
from scipy import stats


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





#---FUNCTIONS TO COMPUTE AN ALPHA COEFFICIENT THAT DEPENDS ON THE POSITION---


def constant_alpha(x, y, value=0.7, **kwargs):
    """ Returns a constant alpha value across the grid. """
    # Ensure alpha is within valid CTRW range (0 < alpha < 1)
    alpha_clipped = np.clip(value, 0.01, 0.99) # Avoid exact 0 or 1
    return np.full(x.shape, alpha_clipped, dtype=np.float32)

def gaussian_alpha(x, y, center_alpha=0.5, edge_alpha=0.9, sigma=0.5, **kwargs):
    """
    Alpha varies radially based on a Gaussian profile.
    Lower alpha (longer waits) near the center.
    Higher alpha (shorter waits) further away.
    """
    distance_sq = x**2 + y**2
    # Gaussian weight (1 at center, decays towards 0)
    weight_center = np.exp(-distance_sq / (2 * sigma**2))
    # Linear interpolation between center and edge alpha based on weight
    alpha = weight_center * center_alpha + (1 - weight_center) * edge_alpha
    # Ensure alpha is within valid CTRW range (0 < alpha < 1)
    return np.clip(alpha, 0.01, 0.99).astype(np.float32)

def linear_gradient_alpha(x, y, min_alpha=0.3, max_alpha=0.9, direction='x', **kwargs):
    """ Alpha varies linearly along a specified direction. """
    if direction == 'x':
        coord = x
        min_coord = np.min(x)
        max_coord = np.max(x)
    elif direction == 'y':
        coord = y
        min_coord = np.min(y)
        max_coord = np.max(y)
    else:
        raise ValueError("Direction must be 'x' or 'y'")

    if max_coord <= min_coord: # Avoid division by zero if grid is flat
        return np.full(x.shape, (min_alpha + max_alpha) / 2, dtype=np.float32)

    # Normalize coordinate to range [0, 1]
    normalized_coord = (coord - min_coord) / (max_coord - min_coord)
    # Interpolate alpha
    alpha = min_alpha + normalized_coord * (max_alpha - min_alpha)
    # Ensure alpha is within valid CTRW range (0 < alpha < 1)
    return np.clip(alpha, 0.01, 0.99).astype(np.float32)



#---ONLY SPATIAL DISORDER FUNCTIONS WITHOUT RESTING PROBABILITY---
def biased_towards_origin(x, y, strength=0.5, **kwargs):
    """
    Increases probability of moving towards the origin (0,0).
    Returns [P+x, P-x, P+y, P-y] summing to 1.
    """
    # Calculate base equal probability
    base_prob = 0.25
    probs = np.full(x.shape + (4,), base_prob, dtype=np.float32)

    # Calculate adjustments towards origin (stronger closer to origin)
    dist_sq = x**2 + y**2
    # Avoid division by zero at origin, add small epsilon
    dist = np.sqrt(dist_sq + 1e-9)

    # Calculate bias factors (negative sign moves towards origin)
    bias_x = -strength * (x / dist) * base_prob # Max bias is strength*base_prob
    bias_y = -strength * (y / dist) * base_prob

    # Apply bias: If x>0, decrease P+x, increase P-x. If x<0, increase P+x, decrease P-x.
    probs[..., IDX_P_X] += bias_x
    probs[..., IDX_M_X] -= bias_x
    probs[..., IDX_P_Y] += bias_y
    probs[..., IDX_M_Y] -= bias_y

    # Clip probabilities to ensure they are valid [0, 1]
    probs = np.clip(probs, 0.0, 1.0)

    # --- Renormalize to ensure sum is exactly 1.0 ---
    prob_sum = np.sum(probs, axis=-1, keepdims=True)
    # Avoid division by zero where sum is zero (shouldn't happen with base_prob=0.25)
    prob_sum[prob_sum < 1e-9] = 1.0
    probs /= prob_sum

    return probs.astype(np.float32)

def vortex_flow(x, y, strength=0.8, **kwargs):
    """
    Creates a swirling probability flow around the origin.
    Returns [P+x, P-x, P+y, P-y] summing to 1.
    """
    base_prob = 0.25
    probs = np.full(x.shape + (4,), base_prob, dtype=np.float32)
    dist_sq = x**2 + y**2 + 1e-9 # Avoid division by zero
    dist = np.sqrt(dist_sq)

    # Tangential direction components (counter-clockwise)
    tangent_x = -y / dist
    tangent_y = x / dist

    # Apply bias based on tangential direction
    # If tangent_x > 0, increase P+x, decrease P-x
    bias_x = strength * tangent_x * base_prob
    # If tangent_y > 0, increase P+y, decrease P-y
    bias_y = strength * tangent_y * base_prob

    probs[..., IDX_P_X] += bias_x
    probs[..., IDX_M_X] -= bias_x
    probs[..., IDX_P_Y] += bias_y
    probs[..., IDX_M_Y] -= bias_y

    probs = np.clip(probs, 0.0, 1.0)
    prob_sum = np.sum(probs, axis=-1, keepdims=True)
    prob_sum[prob_sum < 1e-9] = 1.0
    probs /= prob_sum

    return probs.astype(np.float32)

def gaussian_landscape_per_trial_vectorized(xv, yv, mean_rest_prob=0, std_dev_rest_prob=0.4, max_allowed_rest_prob=0.95):
    """
    Generates a landscape of resting probabilities where each grid cell's
    resting probability is drawn from a Gaussian distribution.
    This landscape is generated per trial due to trial-specific random seeding.

    Args:
        xv (np.ndarray): Meshgrid for x-coordinates.
        yv (np.ndarray): Meshgrid for y-coordinates.
        mean_rest_prob (float): Mean of the Gaussian distribution for resting probability.
        std_dev_rest_prob (float): Standard deviation of the Gaussian distribution.
        max_allowed_rest_prob (float): Maximum allowed resting probability (e.g., 0.95)
                                     to ensure walkers can always move.

    Returns:
        np.ndarray: A (ny, nx, 5) array of probabilities [P+x, P-x, P+y, P-y, Prest].
    """
    ny, nx = xv.shape

    # 1. Draw resting probabilities from a Gaussian distribution for each cell
    # Since np.random.seed() is called per trial, this draw will be unique per trial.
    raw_resting_probs = np.random.normal(loc=mean_rest_prob, scale=std_dev_rest_prob, size=(ny, nx))

    # 2. Clip the resting probabilities to be within [0, max_allowed_rest_prob]
    # This ensures probabilities are valid and that there's always some chance to move.
    clipped_resting_probs = np.clip(raw_resting_probs, 0.0, max_allowed_rest_prob)

    # 3. Calculate the remaining probability available for movement
    probability_for_movement = 1.0 - clipped_resting_probs

    # 4. Distribute the movement probability equally among the four directions
    # Ensure prob_each_direction is not negative (though clipping Prest should prevent this)
    prob_each_direction = np.maximum(0.0, probability_for_movement / 4.0)

    # 5. Construct the output probability array
    # Shape: (ny, nx, 5) for [P+x, P-x, P+y, P-y, Prest]
    output_probs = np.zeros((ny, nx, 5), dtype=np.float32)
    output_probs[..., 0] = prob_each_direction  # P_plus_x
    output_probs[..., 1] = prob_each_direction  # P_minus_x
    output_probs[..., 2] = prob_each_direction  # P_plus_y
    output_probs[..., 3] = prob_each_direction  # P_minus_y
    output_probs[..., 4] = clipped_resting_probs # Prest

    # Optional: Sanity check that probabilities sum to 1.0 for each cell
    # total_cell_probs = np.sum(output_probs, axis=-1)
    # if not np.allclose(total_cell_probs, 1.0):
    #     print("Warning: Probabilities in gaussian_landscape_per_trial do not sum to 1.0 for all cells.")
    #     print(f"Min sum: {np.min(total_cell_probs)}, Max sum: {np.max(total_cell_probs)}")

    return output_probs

def uniform_landscape_per_trial_vectorized(xv, yv, min_rest_prob=0.0, max_rest_prob=1):
    """
    Generates a landscape of resting probabilities where each grid cell's
    resting probability is drawn from a Uniform distribution.
    This landscape is generated per trial due to trial-specific random seeding.

    Args:
        xv (np.ndarray): Meshgrid for x-coordinates.
        yv (np.ndarray): Meshgrid for y-coordinates.
        min_rest_prob (float): Minimum value for the uniform distribution of resting probability.
                               Should be between 0 and 1.
        max_rest_prob (float): Maximum value for the uniform distribution of resting probability.
                               Should be between 0 and 1, and >= min_rest_prob.

    Returns:
        np.ndarray: A (ny, nx, 5) array of probabilities [P+x, P-x, P+y, P-y, Prest].
    """
    ny, nx = xv.shape

    # Validate and clip input probabilities to ensure they are within [0, 1]
    # and min_rest_prob <= max_rest_prob
    # (More robust input validation could be done at the config loading stage)
    actual_min_rest = np.clip(min_rest_prob, 0.0, 1.0)
    actual_max_rest = np.clip(max_rest_prob, 0.0, 1.0)
    if actual_min_rest > actual_max_rest:
        # If min > max after clipping (e.g. min=0.7, max=0.5), swap them or set to a valid range.
        # Here, we'll ensure min is not greater than max.
        actual_min_rest, actual_max_rest = sorted([actual_min_rest, actual_max_rest])
        print(f"Warning: min_rest_prob ({min_rest_prob}) was > max_rest_prob ({max_rest_prob})."
              f" Adjusted to min={actual_min_rest}, max={actual_max_rest}")


    # 1. Draw resting probabilities from a Uniform distribution for each cell.
    # Since np.random.seed() is called per trial, this draw will be unique per trial.
    resting_probs_for_cells = np.random.uniform(low=actual_min_rest, high=actual_max_rest, size=(ny, nx))

    # 2. Calculate the remaining probability available for movement
    probability_for_movement = 1.0 - resting_probs_for_cells

    # 3. Distribute the movement probability equally among the four directions
    # Ensure prob_each_direction is not negative.
    prob_each_direction = np.maximum(0.0, probability_for_movement / 4.0)

    # 4. Construct the output probability array
    # Shape: (ny, nx, 5) for [P+x, P-x, P+y, P-y, Prest]
    output_probs = np.zeros((ny, nx, 5), dtype=np.float32)
    output_probs[..., 0] = prob_each_direction      # P_plus_x
    output_probs[..., 1] = prob_each_direction      # P_minus_x
    output_probs[..., 2] = prob_each_direction      # P_plus_y
    output_probs[..., 3] = prob_each_direction      # P_minus_y
    output_probs[..., 4] = resting_probs_for_cells  # Prest

    return output_probs


def truncated_gaussian_landscape_per_trial_vectorized(xv, yv, mean_rest_prob=0.5, std_dev_rest_prob=0.2, **kwargs):
    """
    Generates a landscape of resting probabilities where each grid cell's
    resting probability is drawn from a Gaussian distribution truncated to [0, 1].
    This landscape is generated per trial due to trial-specific random seeding.

    Args:
        xv (np.ndarray): Meshgrid for x-coordinates.
        yv (np.ndarray): Meshgrid for y-coordinates.
        mean_rest_prob (float): Mean of the underlying Gaussian distribution
                                before truncation.
        std_dev_rest_prob (float): Standard deviation of the underlying Gaussian
                                   distribution before truncation.
        **kwargs: Catches unused parameters.

    Returns:
        np.ndarray: A (ny, nx, 5) array of probabilities [P+x, P-x, P+y, P-y, Prest].
    """
    ny, nx = xv.shape

    # Define the desired clip bounds for the resting probability
    clip_min = 0.0
    clip_max = 1.0

    if std_dev_rest_prob <= 1e-7: # Effectively zero standard deviation
        # If std_dev is very small, truncnorm is ill-defined / just returns the mean.
        # So, we directly use the mean, clipped to the [0, 1] interval.
        resting_probs_for_cells = np.full((ny, nx), np.clip(mean_rest_prob, clip_min, clip_max))
    else:
        # Calculate the 'a' and 'b' parameters for truncnorm.
        # These define the truncation range (clip_min, clip_max) in terms of
        # standard deviations from the mean of the underlying Gaussian.
        # a = (desired_min - mean) / std_dev
        # b = (desired_max - mean) / std_dev
        a_trunc = (clip_min - mean_rest_prob) / std_dev_rest_prob
        b_trunc = (clip_max - mean_rest_prob) / std_dev_rest_prob

        # 1. Draw resting probabilities from a Truncated Gaussian distribution.
        # Since np.random.seed() is called per trial, this draw will be unique.
        # The .rvs method draws random variates.
        resting_probs_for_cells = stats.truncnorm.rvs(
            a_trunc, b_trunc,
            loc=mean_rest_prob,       # Mean of the underlying Gaussian
            scale=std_dev_rest_prob,  # Standard deviation of the underlying Gaussian
            size=(ny, nx)
        )
        # stats.truncnorm.rvs should already return values within [clip_min, clip_max]
        # but an explicit clip can be a safeguard against tiny floating point errors
        # or extreme parameter choices, though it shouldn't strictly be necessary here.
        # resting_probs_for_cells = np.clip(resting_probs_for_cells, clip_min, clip_max)


    # 2. Calculate the remaining probability available for movement
    probability_for_movement = 1.0 - resting_probs_for_cells

    # 3. Distribute the movement probability equally among the four directions
    # Ensure prob_each_direction is not negative.
    prob_each_direction = np.maximum(0.0, probability_for_movement / 4.0) # 4.0 for NUM_DIRECTIONS

    # 4. Construct the output probability array
    # Shape: (ny, nx, 5) for [P+x, P-x, P+y, P-y, Prest]
    output_probs = np.zeros((ny, nx, 5), dtype=np.float32)
    output_probs[..., 0] = prob_each_direction      # P_plus_x
    output_probs[..., 1] = prob_each_direction      # P_minus_x
    output_probs[..., 2] = prob_each_direction      # P_plus_y
    output_probs[..., 3] = prob_each_direction      # P_minus_y
    output_probs[..., 4] = resting_probs_for_cells.astype(np.float32) # Prest, ensure float32

    return output_probs