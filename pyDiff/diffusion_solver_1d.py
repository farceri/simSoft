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
run_static_plots = True  # Set to False to bypass static plots
run_animation_flag = True  # Set to False to bypass animation generation and display
save_animation_gif = False  # Set to True to attempt saving GIF (requires ImageMagick)
save_animation_mp4 = False  # Set to True to attempt saving MP4 (requires FFmpeg)


# --- Numba JIT-compiled function for the core simulation loop ---
@numba.njit(cache=True)
def _run_simulation_loops_numba(p_array, D_x_vals_array, dt_val, dx_val, Nx_val, Nt_val):
    dx2 = dx_val * dx_val
    dt_div_dx2 = dt_val / dx2
    Q_current_buffer = np.empty(Nx_val, dtype=np.float64)

    for j in range(0, Nt_val - 1):
        for k in range(Nx_val):
            Q_current_buffer[k] = D_x_vals_array[k] * p_array[k, j]

        for i in range(1, Nx_val - 1):
            p_array[i, j + 1] = p_array[i, j] + dt_div_dx2 * \
                                (Q_current_buffer[i + 1] - 2 * Q_current_buffer[i] + Q_current_buffer[i - 1])
        p_array[0, j + 1] = p_array[0, j] + dt_div_dx2 * (Q_current_buffer[1] - Q_current_buffer[0])
        p_array[Nx_val - 1, j + 1] = p_array[Nx_val - 1, j] + dt_div_dx2 * \
                                     (Q_current_buffer[Nx_val - 2] - Q_current_buffer[Nx_val - 1])


def solve_pde_ftcs_conservative_numba(D_star, sigma_p, L, T, Nx, Nt):
    x = np.linspace(-L / 2, L / 2, Nx, dtype=np.float64)
    t = np.linspace(0, T, Nt, dtype=np.float64)
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    if np.isclose(sigma_p, 0):
        D_x_vals = D_star * np.ones(Nx, dtype=np.float64)
    else:
        gaussian_term = (1.0 / (np.sqrt(2 * np.pi) * sigma_p)) * \
                        np.exp(-(x ** 2) / (2 * sigma_p ** 2))
        D_x_vals = D_star * (1.0 - gaussian_term)
        D_x_vals = np.maximum(D_x_vals, 1e-9)

    alpha = np.max(D_x_vals) * dt / dx ** 2
    print(f"Stability parameter alpha (based on max(D(x)) dt/dx^2) = {alpha:.4f}")
    if alpha > 0.5:
        print(f"Warning: Stability condition (max(D(x)) dt/dx^2 <= 0.5) may not be met (alpha = {alpha:.4f}).")

    p = np.zeros((Nx, Nt), dtype=np.float64)
    center_index = np.argmin(np.abs(x - 0.0))
    p[center_index, 0] = 1.0 / dx

    _run_simulation_loops_numba(p, D_x_vals, dt, dx, Nx, Nt)
    return x, t, p, D_x_vals


# --- Simulation Parameters ---
D_star_val = 1.0
sigma_p_val = 0.4
L_val = 2.0
T_val = 1.0
Nx_val = 101
Nt_val = 10000  # Using a larger value for a more detailed simulation

# --- Run the simulation ---
print("Running simulation with Numba...")
start_time = timer.time()
x_sol, t_sol, p_sol, Dx_plot = solve_pde_ftcs_conservative_numba(D_star_val, sigma_p_val, L_val, T_val, Nx_val, Nt_val)
end_time = timer.time()
print(f"Simulation completed in {end_time - start_time:.4f} seconds.")

# --- Static Plots (conditionally) ---
if run_static_plots:
    print("Generating static plots...")
    plt.figure(figsize=(12, 12))
    plt.subplot(4, 1, 1)
    plt.plot(x_sol, Dx_plot)
    plt.title(f'Diffusion Coefficient D(x) (D*={D_star_val}, $\sigma_p$={sigma_p_val})')
    plt.xlabel('x');
    plt.ylabel('D(x)');
    plt.grid(True)

    plt.subplot(4, 1, 2)
    time_indices_to_plot = [0, int(Nt_val * 0.1), int(Nt_val * 0.2), int(Nt_val * 0.5), Nt_val - 1]
    for k, time_idx in enumerate(time_indices_to_plot):
        if time_idx < Nt_val:
            plt.plot(x_sol, p_sol[:, time_idx], label=f't = {t_sol[time_idx]:.2f}')
    plt.title('Concentration p(x,t) at different times (Conservative form, Numba)')
    plt.xlabel('x');
    plt.ylabel('p(x,t)');
    plt.legend();
    plt.grid(True)

    ax_3d = plt.subplot(4, 1, 3, projection='3d')
    X_grid, T_grid = np.meshgrid(x_sol, t_sol)
    P_grid = p_sol.T
    stride_t = max(1, Nt_val // 100);
    stride_x = max(1, Nx_val // 50)
    surf = ax_3d.plot_surface(X_grid[::stride_t, ::stride_x], T_grid[::stride_t, ::stride_x],
                              P_grid[::stride_t, ::stride_x], cmap='viridis', edgecolor='none')
    ax_3d.set_title('Surface plot of p(x,t)');
    ax_3d.set_xlabel('x');
    ax_3d.set_ylabel('t');
    ax_3d.set_zlabel('p(x,t)')

    plt.subplot(4, 1, 4)
    integral_p = np.zeros(Nt_val)
    dx_val = x_sol[1] - x_sol[0]
    for j_idx in range(Nt_val):
        integral_p[j_idx] = np.sum(p_sol[:, j_idx]) * dx_val
    plt.plot(t_sol, integral_p)
    plt.title('Integral of p(x,t) over x vs. Time (Conservation Check)')
    plt.xlabel('Time t');
    plt.ylabel('$\int p(x,t) dx$')
    min_integral_val = np.min(integral_p);
    max_integral_val = np.max(integral_p)  # Renamed to avoid conflict
    if min_integral_val > 0 and max_integral_val > 0:
        plt.ylim(min(0.95, 0.99 * min_integral_val), max(1.05, 1.01 * max_integral_val))
    else:
        plt.ylim(0.0, 1.5 if max_integral_val < 1.5 else max_integral_val * 1.1)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"Initial integral of p: {integral_p[0]:.6f}")
    print(f"Final integral of p: {integral_p[-1]:.6f}")
    print(f"Mean integral: {np.mean(integral_p):.6f}, Std dev: {np.std(integral_p):.6e}")
else:
    print("Static plots bypassed.")

# --- Animation Setup (conditionally) ---
if run_animation_flag:
    print("Setting up animation...")
    # Use fewer frames for animation if Nt_val is large for faster generation/display
    animation_display_frames = 200  # Number of frames to display/save in the animation
    if Nt_val < animation_display_frames:  # Ensure we don't try to make more frames than available time steps
        animation_display_frames = Nt_val
    frame_step = max(1, Nt_val // animation_display_frames)

    fig_anim, ax_anim = plt.subplots()
    line, = ax_anim.plot(x_sol, p_sol[:, 0], lw=2)
    ax_anim.set_xlabel('x')
    ax_anim.set_ylabel('p(x,t)')
    ax_anim.set_title('Time evolution of p(x,t) (Numba optimized)')
    ax_anim.grid(True)

    # Dynamic y-axis limit setting
    y_max_anim = 0
    if Nt_val > 1:
        # Find max p value after the initial peak (e.g. after first 1% of steps, or at least 1 step)
        start_idx_for_ymax = min(max(1, Nt_val // 100), Nt_val - 1)
        y_max_anim = np.max(p_sol[:, start_idx_for_ymax:]) * 1.1
    if y_max_anim <= 0:  # Fallback if all values are zero or negative (unlikely)
        y_max_anim = np.max(p_sol[:, 0]) * 0.2 if np.max(p_sol[:, 0]) > 0 else 1.0
    ax_anim.set_ylim(0, y_max_anim)

    time_text = ax_anim.text(0.05, 0.9, '', transform=ax_anim.transAxes)


    def init_animation():
        line.set_ydata(p_sol[:, 0])
        time_text.set_text('')
        return line, time_text


    animation_frame_indices = range(0, Nt_val, frame_step)
    actual_num_animation_frames = len(animation_frame_indices)


    def animate(i_anim_frame):
        actual_time_index = animation_frame_indices[i_anim_frame]
        line.set_ydata(p_sol[:, actual_time_index])
        time_text.set_text(
            f'Time = {t_sol[actual_time_index]:.3f} s (Frame {i_anim_frame + 1}/{actual_num_animation_frames})')
        return line, time_text


    ani = animation.FuncAnimation(fig_anim, animate, frames=actual_num_animation_frames,
                                  init_func=init_animation, blit=True, interval=50, repeat=False)
    plt.show()  # Display the animation

    if save_animation_gif:
        print("Attempting to save GIF...")
        try:
            ani.save('diffusion_animation_numba.gif', writer='pillow',
                     fps=15)  # 'pillow' is a good alternative to 'imagemagick'
            print("GIF saved as diffusion_animation_numba.gif")
        except Exception as e:
            print(f"Could not save GIF: {e}")
            print("Make sure Pillow is installed (`pip install Pillow`) or ImageMagick is available.")

    if save_animation_mp4:
        print("\nAttempting to save MP4...")
        try:
            ani.save('diffusion_animation_numba.mp4', writer='ffmpeg', fps=15, dpi=150)
            print("MP4 saved as diffusion_animation_numba.mp4")
        except Exception as e:
            print(f"Could not save MP4: {e}")
            print("Make sure FFmpeg is installed and in your PATH.")
else:
    print("Animation bypassed.")