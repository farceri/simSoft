import numpy as np


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

