import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numba
import time as timer

# --- Control Flags ---
run_static_plots = True
run_animation_flag = False
save_animation_gif = False  # Set to True to attempt saving GIF (requires Pillow)
save_animation_mp4 = False  # Set to True to attempt saving MP4 (requires FFmpeg)
num_snapshots_to_store = 100  # Fewer snapshots for 2D due to larger data per snapshot


def get_analytical_msd_2d(t_array, D_star, sigma_p, C_perturb_strength):
    """
    Calculates the analytical 2D MSD based on the first-order perturbation theory.
    <r^2(t)> = 4*D*t + C * 4*sigma_p * (sigma_p - sqrt(sigma_p^2 + 2*D*t))
    """
    # Ensure t_array is numpy array for vectorized operations
    t_array_np = np.asarray(t_array)

    # Zeroth-order term
    msd_p0 = 4 * D_star * t_array_np

    # First-order correction term
    # Handle t=0 for the correction part to avoid issues if sigma_p=0 initially, though physically t>0
    # The term sqrt(sigma_p^2 + 2*D_star*t) is well-defined for t>=0
    term_inside_sqrt = sigma_p ** 2 + 2 * D_star * t_array_np

    # Ensure we don't take sqrt of negative if by some numerical error term_inside_sqrt is tiny negative
    # (practically, with t_array_np >= 0 and D_star, sigma_p >=0, this shouldn't happen)
    term_inside_sqrt = np.maximum(term_inside_sqrt, 0)

    correction_factor = sigma_p - np.sqrt(term_inside_sqrt)
    msd_p1_contribution = C_perturb_strength * 4 * sigma_p * correction_factor

    total_msd = msd_p0 + msd_p1_contribution
    return total_msd



# --- Numba JIT-compiled function for a single time step in 2D ---
@numba.njit(cache=True)
def _run_one_step_2d_numba(p_curr_arr, p_next_arr, D_xy_vals_array,
                           dt_val, dx_val, dy_val, Nx_val, Ny_val, Q_buffer):
    """
    Calculates a single time step for the 2D PDE.
    p_curr_arr: current probability distribution (Nx, Ny)
    p_next_arr: array to store the next state (Nx, Ny)
    D_xy_vals_array: diffusion coefficient D(x,y) (Nx, Ny)
    Q_buffer: pre-allocated buffer for D*p (Nx, Ny)
    """
    dx2 = dx_val * dx_val
    dy2 = dy_val * dy_val
    dt_div_dx2 = dt_val / dx2
    dt_div_dy2 = dt_val / dy2

    # Calculate Q_current = D(x,y) * p_current(x,y)
    for i in range(Nx_val):
        for j_idx in range(Ny_val):
            Q_buffer[i, j_idx] = D_xy_vals_array[i, j_idx] * p_curr_arr[i, j_idx]

    # Interior points
    for i in range(1, Nx_val - 1):
        for j_idx in range(1, Ny_val - 1):
            term_xx = Q_buffer[i + 1, j_idx] - 2 * Q_buffer[i, j_idx] + Q_buffer[i - 1, j_idx]
            term_yy = Q_buffer[i, j_idx + 1] - 2 * Q_buffer[i, j_idx] + Q_buffer[i, j_idx - 1]
            p_next_arr[i, j_idx] = p_curr_arr[i, j_idx] + dt_div_dx2 * term_xx + dt_div_dy2 * term_yy

    # Boundary Conditions (Finite Volume Style for zero flux of d(Dp)/dn)
    # Edges (excluding corners)
    for i in range(1, Nx_val - 1):  # Top/Bottom edges
        # Bottom edge (j_idx=0)
        term_xx_b = Q_buffer[i + 1, 0] - 2 * Q_buffer[i, 0] + Q_buffer[i - 1, 0]
        term_yy_b = Q_buffer[i, 1] - Q_buffer[i, 0]  # Flux out of y=0 face is zero
        p_next_arr[i, 0] = p_curr_arr[i, 0] + dt_div_dx2 * term_xx_b + dt_div_dy2 * term_yy_b
        # Top edge (j_idx=Ny_val-1)
        term_xx_t = Q_buffer[i + 1, Ny_val - 1] - 2 * Q_buffer[i, Ny_val - 1] + Q_buffer[i - 1, Ny_val - 1]
        term_yy_t = Q_buffer[i, Ny_val - 2] - Q_buffer[i, Ny_val - 1]  # Flux out of y=Ly face is zero
        p_next_arr[i, Ny_val - 1] = p_curr_arr[i, Ny_val - 1] + dt_div_dx2 * term_xx_t + dt_div_dy2 * term_yy_t

    for j_idx in range(1, Ny_val - 1):  # Left/Right edges
        # Left edge (i=0)
        term_xx_l = Q_buffer[1, j_idx] - Q_buffer[0, j_idx]  # Flux out of x=0 face is zero
        term_yy_l = Q_buffer[0, j_idx + 1] - 2 * Q_buffer[0, j_idx] + Q_buffer[0, j_idx - 1]
        p_next_arr[0, j_idx] = p_curr_arr[0, j_idx] + dt_div_dx2 * term_xx_l + dt_div_dy2 * term_yy_l
        # Right edge (i=Nx_val-1)
        term_xx_r = Q_buffer[Nx_val - 2, j_idx] - Q_buffer[Nx_val - 1, j_idx]  # Flux out of x=Lx face is zero
        term_yy_r = Q_buffer[Nx_val - 1, j_idx + 1] - 2 * Q_buffer[Nx_val - 1, j_idx] + Q_buffer[Nx_val - 1, j_idx - 1]
        p_next_arr[Nx_val - 1, j_idx] = p_curr_arr[Nx_val - 1, j_idx] + dt_div_dx2 * term_xx_r + dt_div_dy2 * term_yy_r

    # Corners
    # Bottom-left (0,0)
    p_next_arr[0, 0] = p_curr_arr[0, 0] + \
                       dt_div_dx2 * (Q_buffer[1, 0] - Q_buffer[0, 0]) + \
                       dt_div_dy2 * (Q_buffer[0, 1] - Q_buffer[0, 0])
    # Top-left (0,Ny-1)
    p_next_arr[0, Ny_val - 1] = p_curr_arr[0, Ny_val - 1] + \
                                dt_div_dx2 * (Q_buffer[1, Ny_val - 1] - Q_buffer[0, Ny_val - 1]) + \
                                dt_div_dy2 * (Q_buffer[0, Ny_val - 2] - Q_buffer[0, Ny_val - 1])
    # Bottom-right (Nx-1,0)
    p_next_arr[Nx_val - 1, 0] = p_curr_arr[Nx_val - 1, 0] + \
                                dt_div_dx2 * (Q_buffer[Nx_val - 2, 0] - Q_buffer[Nx_val - 1, 0]) + \
                                dt_div_dy2 * (Q_buffer[Nx_val - 1, 1] - Q_buffer[Nx_val - 1, 0])
    # Top-right (Nx-1,Ny-1)
    p_next_arr[Nx_val - 1, Ny_val - 1] = p_curr_arr[Nx_val - 1, Ny_val - 1] + \
                                         dt_div_dx2 * (Q_buffer[Nx_val - 2, Ny_val - 1] - Q_buffer[
        Nx_val - 1, Ny_val - 1]) + \
                                         dt_div_dy2 * (Q_buffer[Nx_val - 1, Ny_val - 2] - Q_buffer[
        Nx_val - 1, Ny_val - 1])


def solve_pde_2d_memory_optimized_numba(D_star, sigma_p, dip_strength_factor, Lx, Ly, T, Nx, Ny, Nt, num_snapshots):
    """
    Solves the 2D PDE with memory optimization.
    dip_strength_factor: Factor C for the D(x,y) equation's dip.
    """
    # --- Grids and Parameters ---
    x_coords = np.linspace(-Lx / 2, Lx / 2, Nx, dtype=np.float64)
    y_coords = np.linspace(-Ly / 2, Ly / 2, Ny, dtype=np.float64)
    t_all_steps = np.linspace(0, T, Nt, dtype=np.float64)

    dx = x_coords[1] - x_coords[0] if Nx > 1 else Lx
    dy = y_coords[1] - y_coords[0] if Ny > 1 else Ly
    dt = t_all_steps[1] - t_all_steps[0] if Nt > 1 else T

    # --- Diffusion Coefficient D(x,y) ---
    X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')  # D[i,j] corresponds to x[i], y[j]
    # Using the dip_strength_factor as C from D(x,y) = D_star * (1 - C * exp(-(x^2+y^2)/(2*sigma_p^2)))
    # To match 1D's 1/(sqrt(2pi)sigma_p) for sigma_p=0.5 (~0.798), we can pass this value as dip_strength_factor
    if np.isclose(sigma_p, 0):  # Avoid division by zero if sigma_p is in C
        D_xy_profile = D_star * np.ones((Nx, Ny), dtype=np.float64)
    else:
        exponent = -(X ** 2 + Y ** 2) / (2 * sigma_p ** 2)
        gaussian_part = np.exp(exponent)
        D_xy_profile = D_star * (1.0 - dip_strength_factor * gaussian_part)
        D_xy_profile = np.maximum(D_xy_profile, 1e-9)  # Ensure D is positive

    # --- Stability Check ---
    max_D_val = np.max(D_xy_profile)
    stability_val = 0.0
    if dx > 0 and dy > 0:  # Avoid division by zero
        stability_val = dt * (max_D_val / dx ** 2 + max_D_val / dy ** 2)

    print(f"Parameters: Nx={Nx}, Ny={Ny}, Nt={Nt}, T={T}, Lx={Lx}, Ly={Ly}")
    print(f"dx={dx:.2e}, dy={dy:.2e}, dt={dt:.2e}")
    print(f"Max D(x,y) = {max_D_val:.4f}")
    print(f"Stability value (dt*(maxD/dx^2 + maxD/dy^2)): {stability_val:.4f} (should be <= 0.5)")
    if stability_val > 0.5 and Nt > 1:
        print(f"Warning: Stability condition (<= 0.5) may not be met.")

    # --- Initialization for Simulation Loop ---
    p_current = np.zeros((Nx, Ny), dtype=np.float64)
    p_next = np.zeros((Nx, Ny), dtype=np.float64)
    Q_buffer_main = np.zeros((Nx, Ny), dtype=np.float64)

    # Initial condition: 2D Dirac delta at (0,0)
    center_x_idx = np.argmin(np.abs(x_coords - 0.0))
    center_y_idx = np.argmin(np.abs(y_coords - 0.0))
    if dx > 0 and dy > 0:
        p_current[center_x_idx, center_y_idx] = 1.0 / (dx * dy)
    elif Nx == 1 and Ny == 1:  # Single cell
        p_current[0, 0] = 1.0 / (Lx * Ly)  # Normalize to integral 1

    # --- Setup Snapshot Storage ---
    actual_num_snapshots = min(max(1, num_snapshots), Nt if Nt > 0 else 1)
    if Nt == 0:
        snapshot_time_indices = np.array([0], dtype=int)
    elif actual_num_snapshots == 1:
        snapshot_time_indices = np.array([Nt - 1], dtype=int)  # Store only the last step
    else:
        snapshot_time_indices = np.unique(np.linspace(0, Nt - 1, actual_num_snapshots, dtype=int))
    actual_num_snapshots = len(snapshot_time_indices)

    p_snapshots = np.zeros((Nx, Ny, actual_num_snapshots), dtype=np.float64)  # (Nx, Ny, num_snaps)
    t_snapshots = np.zeros(actual_num_snapshots, dtype=np.float64)
    snapshot_write_ptr = 0
    next_snapshot_idx_to_capture = 0

    if Nt == 0:
        print("Warning: Nt is 0. Storing initial state if requested.")
        if actual_num_snapshots > 0 and snapshot_time_indices[0] == 0:
            p_snapshots[..., 0] = p_current
            t_snapshots[0] = 0.0
        return x_coords, y_coords, t_snapshots, p_snapshots, D_xy_profile

    if next_snapshot_idx_to_capture < actual_num_snapshots and snapshot_time_indices[next_snapshot_idx_to_capture] == 0:
        p_snapshots[..., snapshot_write_ptr] = p_current
        t_snapshots[snapshot_write_ptr] = t_all_steps[0]
        snapshot_write_ptr += 1
        next_snapshot_idx_to_capture += 1

    if Nt > 1:
        print(f"Starting 2D simulation loop for {Nt - 1} time steps...")
        sim_loop_start_time = timer.time()
        for j_time_step in range(0, Nt - 1):
            _run_one_step_2d_numba(p_current, p_next, D_xy_profile, dt, dx, dy, Nx, Ny, Q_buffer_main)
            p_current[:, :] = p_next[:, :]  # Update current state

            current_t_idx_in_all_steps = j_time_step + 1
            if next_snapshot_idx_to_capture < actual_num_snapshots and \
                    current_t_idx_in_all_steps == snapshot_time_indices[next_snapshot_idx_to_capture]:
                p_snapshots[..., snapshot_write_ptr] = p_current
                t_snapshots[snapshot_write_ptr] = t_all_steps[current_t_idx_in_all_steps]
                snapshot_write_ptr += 1
                next_snapshot_idx_to_capture += 1

            if (j_time_step + 1) % (max(1, (Nt - 1) // 100)) == 0 or j_time_step == Nt - 2:
                elapsed_loop_time = timer.time() - sim_loop_start_time
                print(
                    f"  Progress: Step {j_time_step + 1}/{Nt - 1} ({(j_time_step + 1) * 100.0 / (Nt - 1 if Nt > 1 else 1):.1f}%) completed. Loop time: {elapsed_loop_time:.2f}s")
    elif Nt == 1 and actual_num_snapshots == 1 and snapshot_time_indices[
        0] == 0 and snapshot_write_ptr == 0:  # Store IC if Nt=1
        p_snapshots[..., 0] = p_current
        t_snapshots[0] = t_all_steps[0]

    return x_coords, y_coords, t_snapshots, p_snapshots, D_xy_profile


# --- Simulation Parameters ---
D_star_val = 2.5/1000
sigma_p_val = 0.4
# For consistency with 1D dip strength where C = 1/(sqrt(2pi)sigma_p)
# if sigma_p = 0.5, C approx 0.798
#dip_factor_C = 1.0 / (2 * np.pi * sigma_p_val**2) if sigma_p_val > 1e-9 else 0.0
dip_factor_C = 0.39
Lx_val = 2.0
Ly_val = 2.0
T_val = 10  # Reduced T for quicker 2D demo
Nx_val = 101  # Reduced Nx, Ny for quicker 2D demo
Ny_val = 101
# Nt_val must satisfy stability: dt <= 0.25 * dx^2 / D_max (if dx=dy)
# dx = 2.0/50 = 0.04. dx^2 = 0.0016.
# D_max ~ 1.0. dt <= 0.25 * 0.0016 / 1.0 = 0.0004
# If T_val = 0.1, Nt = T_val/dt = 0.1 / 0.0004 = 250.
# Let's choose Nt_val to be slightly above this for safety or to set dt.
# Or set Nt, and dt will be T/Nt. Then check stability.
Nt_val = 100000  # Adjusted for 2D stability and shorter T.
# For user's 10^7, it would be extremely long.


# --- Run the simulation ---
overall_start_time = timer.time()
print(
    f"Preparing to run 2D simulation with Nx={Nx_val}, Ny={Ny_val}, Nt={Nt_val} (will store {num_snapshots_to_store} snapshots)")
x_sol, y_sol, t_sol_snapshots, p_sol_snapshots, D_xy_plot = \
    solve_pde_2d_memory_optimized_numba(D_star_val, sigma_p_val, dip_factor_C, Lx_val, Ly_val, T_val, Nx_val, Ny_val,
                                        Nt_val, num_snapshots_to_store)
overall_end_time = timer.time()
print(f"Total script time (including setup and simulation): {overall_end_time - overall_start_time:.4f} seconds.")

# --- Static Plots (using stored snapshots) ---
if run_static_plots:
    print("Generating static plots for 2D data...")
    # Plot D(x,y)
    plt.figure(figsize=(14, 10))
    plt.subplot(2, 2, 1)
    plt.imshow(D_xy_plot.T, extent=[-Lx_val / 2, Lx_val / 2, -Ly_val / 2, Ly_val / 2], origin='lower', aspect='auto',
               cmap='viridis')
    plt.colorbar(label='D(x,y)')
    plt.title('Diffusion Coefficient D(x,y)')
    plt.xlabel('x');
    plt.ylabel('y')

    # Plot p(x,y,t) at a few snapshot times
    num_snaps_available = p_sol_snapshots.shape[2]
    snap_indices_to_plot = np.linspace(0, num_snaps_available - 1, min(3, num_snaps_available), dtype=int)

    for k_plot, snap_idx in enumerate(snap_indices_to_plot):
        plt.subplot(2, 2, 2 + k_plot)
        plt.imshow(p_sol_snapshots[:, :, snap_idx].T, extent=[-Lx_val / 2, Lx_val / 2, -Ly_val / 2, Ly_val / 2],
                   origin='lower', aspect='auto', cmap='viridis')
        plt.colorbar(label='p(x,y,t)')
        plt.title(f'p(x,y,t) at t={t_sol_snapshots[snap_idx]:.3f}')
        plt.xlabel('x');
        plt.ylabel('y')
    plt.tight_layout()
    plt.show()

    # Conservation Check Plot
    integral_p_snapshots = np.zeros(len(t_sol_snapshots))
    dx_val = x_sol[1] - x_sol[0] if len(x_sol) > 1 else Lx_val
    dy_val = y_sol[1] - y_sol[0] if len(y_sol) > 1 else Ly_val
    for j_idx in range(len(t_sol_snapshots)):
        integral_p_snapshots[j_idx] = np.sum(p_sol_snapshots[:, :, j_idx]) * dx_val * dy_val

    plt.figure(figsize=(7, 5))
    plt.plot(t_sol_snapshots, integral_p_snapshots, marker='o', linestyle='-')
    plt.title('Integral of p(x,y,t) over domain vs. Time (Snapshots)')
    plt.xlabel('Time t (snapshots)');
    plt.ylabel('$\iint p(x,y,t) dx dy$')
    if len(integral_p_snapshots) > 0:
        min_int = np.min(integral_p_snapshots);
        max_int = np.max(integral_p_snapshots)
        if min_int > 0.01 and max_int > 0.01:
            plt.ylim(min(0.95 * min_int, 0.95), max(1.05 * max_int, 1.05))
        else:
            plt.ylim(0.0, max(1.5, max_int * 1.1) if max_int > 0 else 1.5)
    plt.grid(True);
    plt.tight_layout();
    plt.show()
    if len(integral_p_snapshots) > 0:
        print(f"Initial integral (snapshots): {integral_p_snapshots[0]:.6f}")
        print(f"Final integral (snapshots): {integral_p_snapshots[-1]:.6f}")
        print(
            f"Mean integral (snapshots): {np.mean(integral_p_snapshots):.6f}, Std dev: {np.std(integral_p_snapshots):.6e}")

# --- Calculate and Plot 2D MSD ---
if p_sol_snapshots.size > 0:  # Check if there's data
    print("Calculating 2D MSD...")
    dx_val = x_sol[1] - x_sol[0] if len(x_sol) > 1 else Lx_val
    dy_val = y_sol[1] - y_sol[0] if len(y_sol) > 1 else Ly_val

    num_actual_snapshots = p_sol_snapshots.shape[2]
    msd_values = np.zeros(num_actual_snapshots)
    mean_x_sq_values = np.zeros(num_actual_snapshots)  # For <x^2>
    mean_y_sq_values = np.zeros(num_actual_snapshots)  # For <y^2>
    mean_x_vals = np.zeros(num_actual_snapshots)
    mean_y_vals = np.zeros(num_actual_snapshots)

    X_grid, Y_grid = np.meshgrid(x_sol, y_sol, indexing='ij')  # X[i,j]=x[i], Y[i,j]=y[j]
    R_squared_grid = X_grid ** 2 + Y_grid ** 2

    for k in range(num_actual_snapshots):
        current_p_snapshot = p_sol_snapshots[:, :, k]
        msd_values[k] = np.sum(R_squared_grid * current_p_snapshot) * dx_val * dy_val
        mean_x_vals[k] = np.sum(X_grid * current_p_snapshot) * dx_val * dy_val
        mean_y_vals[k] = np.sum(Y_grid * current_p_snapshot) * dx_val * dy_val
        # For verification if needed:
        # mean_x_sq_values[k] = np.sum(X_grid**2 * current_p_snapshot) * dx_val * dy_val
        # mean_y_sq_values[k] = np.sum(Y_grid**2 * current_p_snapshot) * dx_val * dy_val

    print(f"Mean <x> (snapshots): min={np.min(mean_x_vals):.2e}, max={np.max(mean_x_vals):.2e}")
    print(f"Mean <y> (snapshots): min={np.min(mean_y_vals):.2e}, max={np.max(mean_y_vals):.2e}")

    valid_indices_for_plot = []
    if num_actual_snapshots > 0:
        first_valid_idx = 0
        while first_valid_idx < num_actual_snapshots and t_sol_snapshots[first_valid_idx] <= 1e-9:
            first_valid_idx += 1
        if first_valid_idx < num_actual_snapshots:
            valid_indices_for_plot = range(first_valid_idx, num_actual_snapshots)

    t_plot = t_sol_snapshots[valid_indices_for_plot]
    msd_plot = msd_values[valid_indices_for_plot]

    if len(t_plot) > 0:  # Ensure t_plot is not empty
        analytical_msd_values = get_analytical_msd_2d(t_plot, D_star_val, sigma_p_val, dip_factor_C)
        analytical_normalized_msd = analytical_msd_values / t_plot  # t_plot is already filtered for t > 0
    else:
        analytical_msd_values = np.array([])
        analytical_normalized_msd = np.array([])

    if len(t_plot) > 0:
        normalized_msd = msd_plot / t_plot

        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.plot(t_plot, normalized_msd, marker='o', linestyle='-')
        if analytical_normalized_msd.size > 0:
            plt.plot(t_plot, analytical_normalized_msd, marker='x', linestyle='--', color='red',
                     label='Analytical MSD/t (Perturbative)')
        plt.xlabel('Time t');
        plt.ylabel(r'$\langle r^2(t) \rangle / t$');
        plt.title('Normalized MSD vs. Time (2D)')
        plt.grid(True)

        plt.subplot(1, 2, 2)
        positive_msd_log = msd_plot > 1e-9
        if np.any(positive_msd_log):  # t_plot is already positive here
            plt.plot(t_plot[positive_msd_log], msd_plot[positive_msd_log], marker='o', linestyle='-')
            if analytical_msd_values.size > 0:
                # Ensure analytical values are also positive for log plot
                positive_analytical_msd_log = analytical_msd_values > 1e-9
                # Plot only where both t_plot and analytical_msd_values are positive for log scale
                # and align with the t_plot points used for simulation data
                valid_analytical_indices = positive_msd_log & positive_analytical_msd_log
                if np.any(valid_analytical_indices):
                    plt.plot(t_plot[valid_analytical_indices], analytical_msd_values[valid_analytical_indices],
                             marker='x', linestyle='--', color='red', label='Analytical MSD (Perturbative)')

            plt.xscale('log');
            plt.yscale('log')
            plt.xlabel('Time t (log scale)');
            plt.ylabel(r'$\langle r^2(t) \rangle$ (MSD, log scale)')
            plt.title('MSD vs. Time (Log-Log, 2D)');
            plt.grid(True, which="both", ls="-")
        else:
            plt.text(0.5, 0.5, "Not enough positive data for log-log MSD plot", ha='center', va='center')
            plt.title('MSD vs. Time (Log-Log, 2D)')
        plt.tight_layout();
        plt.show()
    else:
        print("Not enough valid time points to plot MSD.")

# --- Animation Setup (conditionally) ---
if run_animation_flag and p_sol_snapshots.size > 0 and p_sol_snapshots.shape[2] > 0:
    print("Setting up 2D animation from stored snapshots...")
    fig_anim, ax_anim = plt.subplots()

    # Initial image setup
    # Use p_sol_snapshots[:, :, 0].T because imshow expects (row, col) ~ (y, x)
    # extent defines coordinates: [xmin, xmax, ymin, ymax]
    img = ax_anim.imshow(p_sol_snapshots[:, :, 0].T,
                         extent=[-Lx_val / 2, Lx_val / 2, -Ly_val / 2, Ly_val / 2],
                         origin='lower', aspect='auto', cmap='viridis', animated=True)
    fig_anim.colorbar(img, ax=ax_anim, label='p(x,y,t)')
    ax_anim.set_xlabel('x');
    ax_anim.set_ylabel('y')
    title_anim = ax_anim.set_title(f'Time evolution of p(x,y,t) at t={t_sol_snapshots[0]:.3f}')

    num_animation_actual_frames = p_sol_snapshots.shape[2]


    def animate_2d(i_frame):
        img.set_array(p_sol_snapshots[:, :, i_frame].T)
        title_anim.set_text(f'Time evolution of p(x,y,t) at t={t_sol_snapshots[i_frame]:.3f}')
        return img, title_anim


    ani = animation.FuncAnimation(fig_anim, animate_2d, frames=num_animation_actual_frames,
                                  blit=True, interval=100, repeat=False)  # blit=True requires animate to return artists
    plt.show()

    if save_animation_gif:
        print("Attempting to save 2D GIF from snapshots...")
        try:
            ani.save('diffusion_2d_anim.gif', writer='pillow', fps=10)
            print("GIF saved as diffusion_2d_anim.gif")
        except Exception as e:
            print(f"Could not save GIF: {e}. Ensure Pillow is installed.")

    if save_animation_mp4:
        print("\nAttempting to save 2D MP4 from snapshots...")
        try:
            ani.save('diffusion_2d_anim.mp4', writer='ffmpeg', fps=10, dpi=150)
            print("MP4 saved as diffusion_2d_anim.mp4")
        except Exception as e:
            print(f"Could not save MP4: {e}. Ensure FFmpeg is installed.")
else:
    print("2D Animation bypassed or no snapshots to animate.")