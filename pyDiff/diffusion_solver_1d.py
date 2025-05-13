import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


'''
Created on May 13 2025, by Luca Sfriso

Questo programma risolve l'equazione di diffusione 1d con un coeffieciente di diffusione gaussiano, condizioni al contorno
naturali e una delta di dirac centrata in x=0 come condizione iniziale

Usa le differenze finite con il metodo FTCS

$\partial_{t}{p(x,t)}=D^{*}(1-frac{1}{\sqrt{2\pi}\sigma_{p}}e^{frac{-x^2}{2 \sigma_{p}^2}})\partial^2_{x}p(x,t)$
'''


def solve_pde_ftcs_conservative(D_star, sigma_p, L, T, Nx, Nt):
    """
    Solves the PDE dp/dt = d^2/dx^2 (D(x)p(x,t)) using the FTCS method.

    Parameters:
    D_star (float): Constant D* in the diffusion coefficient.
    sigma_p (float): Parameter sigma_p in the diffusion coefficient.
    L (float): Length of the spatial domain [-L/2, L/2].
    T (float): Total time for the simulation.
    Nx (int): Number of spatial grid points.
    Nt (int): Number of time steps.

    Returns:
    x (numpy.ndarray): Array of spatial points.
    t (numpy.ndarray): Array of time points.
    p (numpy.ndarray): 2D array of the solution p(x,t).
    D_x_vals (numpy.ndarray): Array of diffusion coefficient values D(x).
    """

    # --- Discretization ---
    x = np.linspace(-L / 2, L / 2, Nx)
    t = np.linspace(0, T, Nt)
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    # --- Diffusion Coefficient D(x) ---
    if np.isclose(sigma_p, 0):
        print("Warning: sigma_p is effectively zero. D(x) might behave unexpectedly.")
        D_x_vals = D_star * np.ones(Nx) # Or handle as a pure constant diffusion D*
    else:
        gaussian_term = (1.0 / (np.sqrt(2 * np.pi) * sigma_p)) * np.exp(-x**2 / (2 * sigma_p**2))
        D_x_vals = D_star * (1.0 - gaussian_term)
        # Ensure D(x) is non-negative, critical for physical meaning and stability
        D_x_vals = np.maximum(D_x_vals, 1e-9) # Prevent D(x) = 0 if it causes issues, though max is better.
                                         # Using a small positive number if D_x can be zero.
                                         # Or ensure D_star * (1-gaussian_term) is always positive by choice of D_star, sigma_p

    # --- Stability Check (Courant-Friedrichs-Lewy condition) ---
    # This condition is primarily for D*p_xx type terms.
    # The equation p_t = (D(x)p)_xx = D_xx p + 2 D_x p_x + D p_xx
    # The D p_xx term is usually dominant for stability.
    alpha = np.max(D_x_vals) * dt / dx**2
    print(f"Stability parameter alpha (based on max(D(x)) dt/dx^2) = {alpha:.4f}")
    if alpha > 0.5:
        print(f"Warning: Stability condition (max(D(x)) dt/dx^2 <= 0.5) may not be met (alpha = {alpha:.4f}).")
        print(f"Consider decreasing dt or increasing dx^2.")
        print(f"Required dt <= {0.5 * dx**2 / np.max(D_x_vals):.2e} for current dx and max(D_x).")

    # --- Initialization ---
    p = np.zeros((Nx, Nt))

    # Initial condition: Dirac delta centered at x=0
    center_index = np.argmin(np.abs(x - 0.0))
    p[center_index, 0] = 1.0 / dx
    # Ensure the rest are zero for t=0
    # This is implicitly handled by np.zeros, but being explicit for p[center_index,0] is key.


    # --- Time Stepping (FTCS for conservative form) ---
    # Pre-calculate Q = D(x)p(x,t) at each step for clarity, or do it inline
    Q = np.zeros(Nx)

    for j in range(0, Nt - 1):  # Time loop
        # Calculate Q_i^j = D_i * p_i^j for all i at current time j
        Q_current = D_x_vals * p[:, j]

        # Interior points
        for i in range(1, Nx - 1):
            p[i, j + 1] = p[i, j] + (dt / dx**2) * \
                          (Q_current[i+1] - 2*Q_current[i] + Q_current[i-1])

        # Boundary Conditions (Zero Flux: J = d(Dp)/dx = 0)
        # Based on finite volume: dp/dt = (J_inner - J_outer)/dx
        # J_outer = 0 at domain boundaries.

        # At i = 0 (left boundary): J_outer (J_{-1/2}) = 0
        # dp0/dt = J_{1/2}/dx = ( (Dp)_1 - (Dp)_0 ) / dx^2
        p[0, j + 1] = p[0, j] + (dt / dx**2) * (Q_current[1] - Q_current[0])

        # At i = Nx-1 (right boundary): J_outer (J_{Nx-1/2}) = 0
        # dp_{Nx-1}/dt = -J_{Nx-3/2}/dx = -( (Dp)_{Nx-1} - (Dp)_{Nx-2} ) / dx^2
        #              = ( (Dp)_{Nx-2} - (Dp)_{Nx-1} ) / dx^2
        p[Nx - 1, j + 1] = p[Nx - 1, j] + (dt / dx**2) * \
                           (Q_current[Nx-2] - Q_current[Nx-1])

    return x, t, p, D_x_vals

# --- Simulation Parameters (using your last successful set for stability) ---
D_star_val = 1.0
sigma_p_val = 0.1
L_val = 2.0
T_val = 1.0
Nx_val = 2001
Nt_val = 10000000 # This Nt should be sufficient for stability given previous alpha ~ 0.22

# --- Run the simulation ---
x_sol, t_sol, p_sol, Dx_plot = solve_pde_ftcs_conservative(D_star_val, sigma_p_val, L_val, T_val, Nx_val, Nt_val)

# --- Plotting D(x) ---
plt.figure(figsize=(12, 12)) # Increased height for the new integral plot

plt.subplot(4, 1, 1) # Changed to 4 rows
plt.plot(x_sol, Dx_plot)
plt.title(f'Diffusion Coefficient D(x) (D*={D_star_val}, $\sigma_p$={sigma_p_val})')
plt.xlabel('x')
plt.ylabel('D(x)')
plt.grid(True)

# --- Plotting the solution p(x,t) ---
plt.subplot(4, 1, 2) # Changed to 4 rows
time_indices_to_plot = [0, int(Nt_val / 10), int(Nt_val / 5), int(Nt_val / 2), Nt_val - 1]
for k, time_idx in enumerate(time_indices_to_plot):
    # Check if time_idx is within bounds
    if time_idx < Nt_val:
        plt.plot(x_sol, p_sol[:, time_idx], label=f't = {t_sol[time_idx]:.2f}')
plt.title('Concentration p(x,t) at different times (Conservative form)')
plt.xlabel('x')
plt.ylabel('p(x,t)')
plt.legend()
plt.grid(True)

# --- Surface plot of p(x,t) ---
ax_3d = plt.subplot(4, 1, 3, projection='3d') # Changed to 4 rows
X_grid, T_grid = np.meshgrid(x_sol, t_sol)
P_grid = p_sol.T
stride_t = max(1, Nt_val // 50)
stride_x = max(1, Nx_val // 50)
surf = ax_3d.plot_surface(X_grid[::stride_t, ::stride_x], T_grid[::stride_t, ::stride_x], P_grid[::stride_t, ::stride_x], cmap='viridis', edgecolor='none')
ax_3d.set_title('Surface plot of p(x,t)')
ax_3d.set_xlabel('x')
ax_3d.set_ylabel('t')
ax_3d.set_zlabel('p(x,t)')


# --- Check conservation of total probability (integral of p(x,t) over x) ---
plt.subplot(4, 1, 4) # Changed to 4 rows
integral_p = np.zeros(Nt_val)
dx_val = x_sol[1] - x_sol[0]
for j in range(Nt_val):
    integral_p[j] = np.sum(p_sol[:, j]) * dx_val

plt.plot(t_sol, integral_p)
plt.title('Integral of p(x,t) over x vs. Time (Conservation Check)')
plt.xlabel('Time t')
plt.ylabel('$\int p(x,t) dx$')
# Set y-limits around 1.0 for better visualization of conservation
min_integral = np.min(integral_p)
max_integral = np.max(integral_p)
if min_integral > 0 and max_integral > 0 : # Check if integral values are valid
    plt.ylim(min(0.9, 0.98 * min_integral), max(1.1, 1.02 * max_integral))
else: # Fallback if integrals are zero or negative (should not happen with proper IC)
    plt.ylim(0.0, 1.5 if max_integral < 1.5 else max_integral * 1.1)

plt.grid(True)

plt.tight_layout()
plt.show()

print(f"Initial integral of p: {integral_p[0]:.6f}")
print(f"Final integral of p: {integral_p[-1]:.6f}")
print(f"Mean integral: {np.mean(integral_p):.6f}, Std dev: {np.std(integral_p):.6e}")
