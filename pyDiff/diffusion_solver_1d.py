import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import  numba
import time as timer

'''
Created on May 13 2025, by Luca Sfriso

Questo programma risolve l'equazione di diffusione 1d con un coeffieciente di diffusione gaussiano, condizioni al contorno
naturali e una delta di dirac centrata in x=0 come condizione iniziale

Usa le differenze finite con il metodo FTCS

$\partial_{t}{p(x,t)}=D^{*}(1-frac{1}{\sqrt{2\pi}\sigma_{p}}e^{frac{-x^2}{2 \sigma_{p}^2}})\partial^2_{x}p(x,t)$
'''

# --- Control Flags ---
run_static_plots = True
run_animation_flag = True
save_animation_gif = False
save_animation_mp4 = False
num_snapshots_to_store = 1000


# --- Numba JIT-compiled function for a single time step ---
@numba.njit(cache=True)
def _run_one_step_numba(p_curr_arr, p_next_arr, D_x_vals_array,
                        dt_val, dx_val, Nx_val, Q_buffer):
    dx2 = dx_val * dx_val
    dt_div_dx2 = dt_val / dx2
    for k in range(Nx_val):
        Q_buffer[k] = D_x_vals_array[k] * p_curr_arr[k]
    for i in range(1, Nx_val - 1):
        p_next_arr[i] = p_curr_arr[i] + dt_div_dx2 * \
                        (Q_buffer[i + 1] - 2 * Q_buffer[i] + Q_buffer[i - 1])
    p_next_arr[0] = p_curr_arr[0] + dt_div_dx2 * (Q_buffer[1] - Q_buffer[0])
    p_next_arr[Nx_val - 1] = p_curr_arr[Nx_val - 1] + dt_div_dx2 * \
                             (Q_buffer[Nx_val - 2] - Q_buffer[Nx_val - 1])


def solve_pde_memory_optimized_numba(D_star, sigma_p, L, T, Nx, Nt, num_snapshots):
    x_coords = np.linspace(-L / 2, L / 2, Nx, dtype=np.float64)
    t_all_steps = np.linspace(0, T, Nt, dtype=np.float64)
    dx = x_coords[1] - x_coords[0] if Nx > 1 else L
    dt = t_all_steps[1] - t_all_steps[0] if Nt > 1 else T

    if np.isclose(sigma_p, 0):
        D_x_profile = D_star * np.ones(Nx, dtype=np.float64)
    else:
        gaussian_term = (1.0 / (np.sqrt(2 * np.pi) * sigma_p)) * \
                        np.exp(-(x_coords ** 2) / (2 * sigma_p ** 2))
        D_x_profile = D_star * (1.0 - gaussian_term)
        D_x_profile = np.maximum(D_x_profile, 1e-9)

    alpha = np.max(D_x_profile) * dt / dx ** 2 if dx > 0 else np.inf
    print(f"Parameters: Nx={Nx}, Nt={Nt}, T={T}, L={L}, dx={dx:.2e}, dt={dt:.2e}")
    print(f"Max D(x) = {np.max(D_x_profile):.4f}")
    print(f"Stability parameter alpha (max(D(x)) dt/dx^2) = {alpha:.4f}")
    if alpha > 0.5 and Nt > 1:
        print(f"Warning: Stability condition (alpha <= 0.5) may not be met.")

    p_current = np.zeros(Nx, dtype=np.float64)
    p_next = np.zeros(Nx, dtype=np.float64)
    Q_buffer_main = np.zeros(Nx, dtype=np.float64)

    center_index = np.argmin(np.abs(x_coords - 0.0))
    if dx > 0:
        p_current[center_index] = 1.0 / dx
    elif Nx == 1:
        p_current[center_index] = 1.0

    actual_num_snapshots = min(max(1, num_snapshots), Nt if Nt > 0 else 1)  # Ensure at least 1 if Nt=0 (for IC)
    if Nt == 0:  # Handle Nt=0 specifically for snapshot indices
        snapshot_time_indices = np.array([0], dtype=int)
    elif actual_num_snapshots == 1:
        snapshot_time_indices = np.array([Nt - 1], dtype=int)
    else:
        snapshot_time_indices = np.unique(np.linspace(0, Nt - 1, actual_num_snapshots, dtype=int))

    actual_num_snapshots = len(snapshot_time_indices)

    p_snapshots = np.zeros((Nx, actual_num_snapshots), dtype=np.float64)
    t_snapshots = np.zeros(actual_num_snapshots, dtype=np.float64)
    snapshot_write_ptr = 0
    next_snapshot_idx_to_capture = 0

    if Nt == 0:
        print("Warning: Nt is 0. Storing initial state as the only snapshot.")
        if actual_num_snapshots > 0 and snapshot_time_indices[0] == 0:
            p_snapshots[:, 0] = p_current
            t_snapshots[0] = 0.0
        return x_coords, t_snapshots, p_snapshots, D_x_profile

    if next_snapshot_idx_to_capture < actual_num_snapshots and snapshot_time_indices[next_snapshot_idx_to_capture] == 0:
        p_snapshots[:, snapshot_write_ptr] = p_current
        t_snapshots[snapshot_write_ptr] = t_all_steps[0]
        snapshot_write_ptr += 1
        next_snapshot_idx_to_capture += 1

    if Nt > 1:
        print(f"Starting simulation loop for {Nt - 1} time steps...")
        sim_loop_start_time = timer.time()
        for j in range(0, Nt - 1):
            _run_one_step_numba(p_current, p_next, D_x_profile, dt, dx, Nx, Q_buffer_main)
            p_current[:] = p_next[:]
            current_time_step_index = j + 1
            if next_snapshot_idx_to_capture < actual_num_snapshots and \
                    current_time_step_index == snapshot_time_indices[next_snapshot_idx_to_capture]:
                p_snapshots[:, snapshot_write_ptr] = p_current
                t_snapshots[snapshot_write_ptr] = t_all_steps[current_time_step_index]
                snapshot_write_ptr += 1
                next_snapshot_idx_to_capture += 1
            if (j + 1) % (max(1, (Nt - 1) // 100)) == 0 or j == Nt - 2:
                elapsed_loop_time = timer.time() - sim_loop_start_time
                print(
                    f"  Progress: Step {j + 1}/{Nt - 1} ({(j + 1) * 100.0 / (Nt - 1 if Nt > 1 else 1):.1f}%) completed. Loop time: {elapsed_loop_time:.2f}s")
    elif Nt == 1 and actual_num_snapshots == 1 and snapshot_time_indices[0] == 0 and snapshot_write_ptr == 0:
        p_snapshots[:, 0] = p_current
        t_snapshots[0] = t_all_steps[0]

    return x_coords, t_snapshots, p_snapshots, D_x_profile


# --- Simulation Parameters ---
D_star_val = 1.0
sigma_p_val = 0.4
L_val = 2.0
T_val = 1.0
Nx_val = 101
Nt_val = 10000

# --- Run the simulation ---
overall_start_time = timer.time()
print(f"Preparing to run simulation with Nx={Nx_val}, Nt={Nt_val} (will store {num_snapshots_to_store} snapshots)")
x_sol, t_sol_snapshots, p_sol_snapshots, Dx_plot = \
    solve_pde_memory_optimized_numba(D_star_val, sigma_p_val, L_val, T_val, Nx_val, Nt_val, num_snapshots_to_store)
overall_end_time = timer.time()
print(f"Total script time (including setup and simulation): {overall_end_time - overall_start_time:.4f} seconds.")

# --- Static Plots (conditionally) ---
if run_static_plots:
    print("Generating static plots from stored snapshots...")
    plt.figure(figsize=(12, 12))
    plt.subplot(4, 1, 1)
    plt.plot(x_sol, Dx_plot)
    plt.title(f'Diffusion Coefficient D(x) (D*={D_star_val}, $\sigma_p$={sigma_p_val})')
    plt.xlabel('x');
    plt.ylabel('D(x)');
    plt.grid(True)

    plt.subplot(4, 1, 2)
    num_snapshots_plotted = len(t_sol_snapshots)
    plot_indices = np.linspace(0, num_snapshots_plotted - 1, min(5, num_snapshots_plotted), dtype=int)
    for k_idx in plot_indices:
        plt.plot(x_sol, p_sol_snapshots[:, k_idx], label=f't = {t_sol_snapshots[k_idx]:.3f} (snapshot {k_idx})')
    plt.title(f'Concentration p(x,t) at different snapshot times')
    plt.xlabel('x');
    plt.ylabel('p(x,t)');
    plt.legend();
    plt.grid(True)

    ax_3d = plt.subplot(4, 1, 3, projection='3d')
    if p_sol_snapshots.shape[0] > 1 and p_sol_snapshots.shape[1] > 1:
        X_grid, T_grid_snapshots = np.meshgrid(x_sol, t_sol_snapshots)
        P_grid_snapshots = p_sol_snapshots.T
        # Ensure rcount and ccount are at least 1
        rcount_val = max(1, P_grid_snapshots.shape[0] // (max(1, P_grid_snapshots.shape[0] // 50)))
        ccount_val = max(1, P_grid_snapshots.shape[1] // (max(1, P_grid_snapshots.shape[1] // 50)))
        if P_grid_snapshots.size > 0:
            surf = ax_3d.plot_surface(X_grid, T_grid_snapshots, P_grid_snapshots, cmap='viridis', edgecolor='none',
                                      rcount=rcount_val, ccount=ccount_val)
            ax_3d.set_title('Surface plot of p(x,t) (from snapshots)')
            ax_3d.set_xlabel('x');
            ax_3d.set_ylabel('t');
            ax_3d.set_zlabel('p(x,t)')
        else:
            ax_3d.set_title('Surface plot not shown (strided data empty)')
    else:
        ax_3d.set_title('Surface plot not shown (too few snapshot dimensions)')

    plt.subplot(4, 1, 4)
    integral_p_snapshots = np.zeros(len(t_sol_snapshots))
    dx_val = x_sol[1] - x_sol[0] if len(x_sol) > 1 else 1.0
    for j_idx in range(len(t_sol_snapshots)):
        integral_p_snapshots[j_idx] = np.sum(p_sol_snapshots[:, j_idx]) * dx_val
    plt.plot(t_sol_snapshots, integral_p_snapshots, marker='o', linestyle='-')
    plt.title('Integral of p(x,t) over x vs. Time (from stored snapshots)')
    plt.xlabel('Time t (snapshots)');
    plt.ylabel('$\int p(x,t) dx$')
    if len(integral_p_snapshots) > 0:
        min_integral_val = np.min(integral_p_snapshots);
        max_integral_val = np.max(integral_p_snapshots)
        if min_integral_val > 0.01 and max_integral_val > 0.01:
            plt.ylim(min(0.95 * min_integral_val, 0.95), max(1.05 * max_integral_val, 1.05))
        else:
            plt.ylim(0.0, max(1.5, max_integral_val * 1.1) if max_integral_val > 0 else 1.5)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    if len(integral_p_snapshots) > 0:
        print(f"Initial integral (from snapshots): {integral_p_snapshots[0]:.6f}")
        print(f"Final integral (from snapshots): {integral_p_snapshots[-1]:.6f}")
        print(
            f"Mean integral (snapshots): {np.mean(integral_p_snapshots):.6f}, Std dev: {np.std(integral_p_snapshots):.6e}")
else:
    print("Static plots bypassed.")

# --- Calculate and Plot Second Moment (MSD) ---
print("Calculating second moments (MSD)...")
dx_val = x_sol[1] - x_sol[0] if len(x_sol) > 1 else L_val
x_squared = x_sol ** 2
num_actual_snapshots = p_sol_snapshots.shape[1]
msd_values = np.zeros(num_actual_snapshots)  # Renamed from second_moments for clarity
mean_x_values = np.zeros(num_actual_snapshots)

for k in range(num_actual_snapshots):
    msd_values[k] = np.sum(x_squared * p_sol_snapshots[:, k]) * dx_val
    mean_x_values[k] = np.sum(x_sol * p_sol_snapshots[:, k]) * dx_val

print(
    f"Mean <x> values (should be close to 0): min={np.min(mean_x_values):.2e}, max={np.max(mean_x_values):.2e}, avg={np.mean(mean_x_values):.2e}")

# Prepare data for plots, excluding t=0 or very small t where MSD might be 0 or t is too small for log.
valid_indices_for_plot = []
if num_actual_snapshots > 0:
    first_valid_idx = 0
    # For MSD vs t, we can often include t=0 if MSD(0)=0.
    # But for log-log, t must be > 0 and MSD > 0.
    # We'll filter for t > small_epsilon for log-log t-axis
    # And msd > small_epsilon for log-log msd-axis
    while first_valid_idx < num_actual_snapshots and t_sol_snapshots[
        first_valid_idx] <= 1e-9:  # Ensure time is positive for log scale later
        first_valid_idx += 1

    if first_valid_idx < num_actual_snapshots:
        valid_indices_for_plot = range(first_valid_idx, num_actual_snapshots)

t_plot = t_sol_snapshots[valid_indices_for_plot]
msd_plot = msd_values[valid_indices_for_plot]

if len(t_plot) > 0:
    # Plot 1: MSD/t vs t (linear scale) - Kept this as it might still be useful
    normalized_msd = msd_plot / t_plot  # Recalculate for the valid t_plot

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.plot(t_plot, normalized_msd, marker='o', linestyle='-')
    plt.xlabel('Time t')
    plt.ylabel(r'$\langle x^2(t) \rangle / t$')
    plt.title('Normalized MSD vs. Time')
    plt.grid(True)

    # Plot 2: MSD vs t (log-log scale) - User's request
    plt.subplot(1, 2, 2)
    # Ensure data is positive for log scale for both t_plot and msd_plot
    positive_t_log = t_plot > 1e-9  # Already handled by valid_indices_for_plot selection
    positive_msd_log = msd_plot > 1e-9  # MSD(t=0) is 0, so msd_plot[0] (if t_plot starts at t>0) might be >0

    # Combine conditions for data points to be plotted on log-log
    log_plot_indices = positive_msd_log  # t_plot is already positive here

    if np.any(log_plot_indices):
        plt.plot(t_plot[log_plot_indices],
                 msd_plot[log_plot_indices],
                 marker='o', linestyle='-')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Time t (log scale)')
        plt.ylabel(r'$\langle x^2(t) \rangle$ (MSD, log scale)')  # Updated Y-label
        plt.title('MSD vs. Time (Log-Log)')  # Updated Title
        plt.grid(True, which="both", ls="-")
    else:
        plt.text(0.5, 0.5, "Not enough positive data for log-log MSD plot", ha='center', va='center')
        plt.xlabel('Time t')
        plt.ylabel(r'$\langle x^2(t) \rangle$ (MSD)')
        plt.title('MSD vs. Time (Log-Log)')

    plt.tight_layout()
    plt.show()
else:
    print("Not enough valid time points to plot MSD.")

# --- Animation Setup (conditionally) ---
if run_animation_flag and p_sol_snapshots.shape[1] > 0:
    print("Setting up animation from stored snapshots...")
    fig_anim, ax_anim = plt.subplots()
    line, = ax_anim.plot(x_sol, p_sol_snapshots[:, 0], lw=2)
    ax_anim.set_xlabel('x');
    ax_anim.set_ylabel('p(x,t)')
    ax_anim.set_title('Time evolution of p(x,t) (Memory Optimized, Numba)')
    ax_anim.grid(True)
    y_max_anim = 0
    if p_sol_snapshots.shape[1] > 1:
        start_idx_for_ymax = min(1, p_sol_snapshots.shape[1] - 1)
        y_max_anim = np.max(p_sol_snapshots[:, start_idx_for_ymax:]) * 1.1
    if y_max_anim <= 0: y_max_anim = np.max(p_sol_snapshots[:, 0]) * 0.2 if np.max(p_sol_snapshots[:, 0]) > 0 else 1.0
    ax_anim.set_ylim(0, y_max_anim)
    time_text = ax_anim.text(0.05, 0.9, '', transform=ax_anim.transAxes)
    num_animation_actual_frames = p_sol_snapshots.shape[1]


    def init_animation():
        line.set_ydata(p_sol_snapshots[:, 0])
        time_text.set_text(f'Time = {t_sol_snapshots[0]:.3f} s (Frame 1/{num_animation_actual_frames})')
        return line, time_text


    def animate(i_frame):
        line.set_ydata(p_sol_snapshots[:, i_frame])
        time_text.set_text(
            f'Time = {t_sol_snapshots[i_frame]:.3f} s (Frame {i_frame + 1}/{num_animation_actual_frames})')
        return line, time_text


    ani = animation.FuncAnimation(fig_anim, animate, frames=num_animation_actual_frames,
                                  init_func=init_animation, blit=True, interval=100, repeat=False)
    plt.show()
    if save_animation_gif:
        print("Attempting to save GIF from snapshots...")
        try:
            ani.save('diffusion_mem_opt_anim.gif', writer='pillow', fps=10)
            print("GIF saved as diffusion_mem_opt_anim.gif")
        except Exception as e:
            print(f"Could not save GIF: {e}. Ensure Pillow is installed.")
    if save_animation_mp4:
        print("\nAttempting to save MP4 from snapshots...")
        try:
            ani.save('diffusion_mem_opt_anim.mp4', writer='ffmpeg', fps=10, dpi=150)
            print("MP4 saved as diffusion_mem_opt_anim.mp4")
        except Exception as e:
            print(f"Could not save MP4: {e}. Ensure FFmpeg is installed.")
else:
    print("Animation bypassed or no snapshots to animate.")