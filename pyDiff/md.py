import numpy as np
import sys
import os
from matplotlib import pyplot as plt

class MolecularDynamics:
    def __init__(self, num_particles=100, temperature=0.1, dt=0.001, gamma=1.0, mass=1.0, Lx=1.0, Ly=1.0, interaction=False):
        self.num_particles = num_particles
        self.mass = mass
        self.dt = dt
        self.gamma = gamma  # Friction coefficient for Langevin dynamics
        self.temperature = temperature
        self.kB = 1.0  # Boltzmann constant (arbitrary units)
        self.mean_vel = np.sqrt(self.kB * self.temperature / self.mass)
        self.box_size = np.array([Lx, Ly])
        print(f"Created md object with settings:")
        print(f"Number of particles: {self.num_particles:d}\nTemperature: {self.temperature:.4f}")
        print(f"Time step: {self.dt:.4f}\nBox size: Lx {Lx:.4f} and Ly {Ly:.4f}")
        
        # Initialize positions and velocities randomly
        self.positions = (np.random.rand(num_particles, 2) - 0.5) * self.box_size # Flat distribution in [0,1] and shifted
        self.initial_positions = np.copy(self.positions) # Store initial positions to compute the MSD
        self.velocities = np.random.randn(num_particles, 2) * self.mean_vel # Normal distribution centered in 0 and variance 1
        self.forces = np.zeros((num_particles, 2))

        # Set size for interacting particles - no force is implemented yet
        if interaction:
            """Set particle sizes based on the shortest interparticle distance."""
            min_distance = np.inf
            for i in range(self.num_particles):
                for j in range(i + 1, self.num_particles):
                    r_ij = np.linalg.norm(self.positions[i] - self.positions[j])
                    if r_ij < min_distance:
                        min_distance = r_ij
            # Set all particles' radii to the minimum interparticle distance
            self.radii = np.full(self.num_particles, min_distance)
            print(f'Average particle diameter: {np.mean(self.radii).df}')

    def apply_pbc(self):
        """Apply periodic boundary conditions to keep particles inside the simulation box."""
        self.positions = (self.positions + self.box_size / 2) % self.box_size - self.box_size / 2

    def compute_forces(self):
        """TODO: add forces using Lennard-Jones potential."""
        self.forces = np.zeros((self.num_particles, 2))  # Reset forces
        # (LJ forces could be added here, for now, assume no interparticle forces)

    def langevin_force(self):
        """Compute stochastic white noise and friction forces."""
        noise = np.sqrt(2 * self.kB * self.temperature * self.gamma / self.dt) * np.random.randn(self.num_particles, 2)
        return -self.gamma * self.velocities + noise

    def velocity_verlet_nve(self):
        """Velocity Verlet integration for NVE dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        #self.apply_pbc()
        self.compute_forces()
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def velocity_verlet_langevin(self):
        """Velocity Verlet integration for Langevin dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        #self.apply_pbc()
        self.compute_forces()
        self.forces += self.langevin_force()
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def compute_temperature(self):
        """Compute the temperature of the system from the kinetic energy."""
        kinetic_energy = 0.5 * self.mass * np.sum(self.velocities ** 2)
        return kinetic_energy / self.num_particles
    
    def compute_msd(self):
        """Compute the mean squared displacement."""
        displacement = self.positions - self.initial_positions
        msd = np.mean(np.sum(displacement ** 2, axis=1))
        return msd

# Example usage - line below allows for exporting the md.py file as a package in another script
# Ex.: import md
if __name__ == '__main__':
    # Read input parameters
    directory = sys.argv[1] # Directory for input and output
    integrator = sys.argv[2] # Integrator type - options are NVE and Langevin
    num_particles = int(sys.argv[3])
    temperature = float(sys.argv[4])
    num_steps = int(float(sys.argv[5])) # Number of integration steps
    save_freq = int(num_steps/100)
    print_freq = int(num_steps/10)
    
    # Create md object with input settings - more settings can be added
    md = MolecularDynamics(num_particles, temperature)

    # Create arrays for storing energy and msd
    temp = np.empty(0)
    msd = np.empty(0)
    # Run integration, store and print data at given frequency
    for step in range(num_steps + save_freq):
        if integrator == 'nve':
            md.velocity_verlet_nve()
        elif integrator == 'langevin':
            md.velocity_verlet_langevin()
        if step % save_freq == 0:
            temp = np.append(temp, md.compute_temperature())
            msd = np.append(msd, md.compute_msd())
        if step % print_freq == 0:    
            print(f"Step {step}: Temperature = {temp[-1]:.4f}, MSD = {msd[-1]:.4f}")
    
    # Store time, temperature and msd in a single file
    time = np.arange(0, num_steps + save_freq, save_freq) * md.dt # Define time array
    np.savetxt(directory + os.sep + 'md_data.dat', np.column_stack((time, temp, msd)))

    # Plot energy and msd versus time
    fig, ax = plt.subplots(2, 1, figsize = (7, 7), sharex = True, dpi = 120)
    ax[0].plot(time, msd, color='k', linestyle='solid', marker='o', markersize='4', fillstyle='none')
    ax[1].plot(time, temp, color='k', linewidth=0.9, linestyle='solid', marker='o', markersize='6', fillstyle='none')
    ax[0].tick_params(axis='both', labelsize=14)
    ax[1].tick_params(axis='both', labelsize=14)
    ax[1].set_xlabel("$Simulation$ $time,$ $t$", fontsize=16)
    ax[0].set_ylabel("$MSD$", fontsize=16)
    ax[1].set_ylabel("$Temperature,$ $T$", fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    #plt.savefig("/home/francesco/Pictures/test.png", transparent=False, format = "png")
    plt.show()
