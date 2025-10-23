import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import time
np.random.seed(0)
# python md.py '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test' 'nve' 10 0.1 10000

class MolecularDynamics:

    def __init__(self, num_particles=100, temperature=0.1, dt=0.001, gamma=1.0, mass=1.0, Lx=10, Ly=10, interaction=False):

        self.num_particles = num_particles
        self.mass = mass
        self.dt = dt
        self.gamma = gamma  # Friction coefficient for Langevin dynamics
        self.temperature = temperature
        self.kB = 1.0  # Boltzmann constant (arbitrary units)
        self.box_size = np.array([Lx, Ly])
        self.neighbours = [] # Initializing neighbour list
        self.cutoff = 3 # Cutoff for disk neighbours
        self.division = 10  # How much to devide the total space for cell neighbours
        self.maxDist = np.minimum(self.cutoff, np.sqrt(2*(self.box_size[0]/self.division)**2)) # Necessary for the force shift
        self.sigma = 1  # Sigma value for LJ potential
        self.epsilon = 1  # Epsilon value for LJ potential
        print(f"Created md object with settings:")
        print(f"Number of particles: {self.num_particles:d}\nTemperature: {self.temperature:.4f}")
        print(f"Time step: {self.dt:.4f}\nBox size: Lx {Lx:.4f} and Ly {Ly:.4f}")
        
        # Initialize positions in a grid 
        nSide = int(np.ceil(np.sqrt(self.num_particles)))
        # Create grid positions
        x = (np.arange(nSide) + 0.5) * (self.box_size[0]/nSide) - (self.box_size[0]/nSide)/2
        y = (np.arange(nSide) + 0.5) * (self.box_size[1]/nSide) - (self.box_size[1]/nSide)/2
        xv, yv = np.meshgrid(x, y)
        positions = np.vstack([xv.ravel(), yv.ravel()]).T
        # Only return the first num_particles if grid has extra points
        self.positions = positions[:self.num_particles]

        # Initialize positions
        self.initial_positions = np.copy(self.positions) # Store initial positions to compute the MSD
        self.velocities = np.random.normal(0, 1, (self.num_particles, 2)) # Maxwell-Boltmann
        self.velocities = self.velocities - (np.sum(self.velocities, axis=0)/self.num_particles) # Remove center of mass
        # Re-sample outliers
        vMax = (self.maxDist/2)
        outliers = np.sqrt(np.sum(self.velocities**2, axis=1)) > vMax
        while outliers.any():
            self.velocities[outliers] = np.random.normal(0, 1, (np.sum(outliers), 2))
            self.velocities = self.velocities - (np.sum(self.velocities, axis=0)/self.num_particles) # Remove center of mass
            outliers = np.sqrt(np.sum(self.velocities**2, axis=1)) > vMax
        self.velocities = self.velocities * np.sqrt((self.num_particles * self.temperature)/(0.5 * self.mass * np.sum(self.velocities ** 2))) # Scale to have initial temperature
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

    def compute_cell_neighbours(self):
        """"Computing nearest neighbours based on cell subdivision."""

        self.maxDist = np.sqrt(2*(self.box_size[0]/self.division)**2)
        self.neighbours = [] # Reset neighbours
        head = np.zeros((self.division, self.division), dtype=int)
        cell = np.zeros((self.num_particles, 2), dtype=int)
        list = np.zeros(self.num_particles, dtype=int)
        near_cells = [(1, 1), (1, 0), (1, -1), (0, 1)]
        other_cell = [0, 0]

        for ii in range(self.num_particles):
            cell[ii, :] = np.floor((self.positions[ii, :] + 0.5)*(self.division)) # which cell ii belongs to
            list[ii] = head[cell[ii, 0], cell[ii, 1]] # point ii to previous head of the cell, 0 if it is first in cell
            head[cell[ii, 0], cell[ii, 1]] = ii # ii is now the new head of the cell

        for ii in range(num_particles):
            ii_list = []
            current_cell = (cell[ii, 0], cell[ii, 1])
            current_head = head[cell[ii, 0], cell[ii, 1]]

            # Neighbours in current cell
            while current_head != 0: # if there are other neighbours...
                if current_head != ii: # ...and it's not yourself
                    ii_list.append(current_head) # append neighbour
                current_head = list[current_head] # next head in line

            # Neighbours in near cells
            for value in near_cells:
                other_cell[0] = current_cell[0] + value[0]
                other_cell[1] = current_cell[1] + value[1]
                current_head = head[other_cell[0]%self.division, other_cell[1]%self.division]
                while current_head != 0: # if there are other neighbours...
                    ii_list.append(current_head) # append neighbour
                    current_head = list[current_head] # next head in line
                    
            self.neighbours.append(ii_list.copy())

    def compute_disk_neighbours(self):
        """"Computing nearest neighbours based on simple disk distance."""

        self.neighbours = [] # Reset neighbours
        self.maxDist = self.cutoff

        for ii in range(self.num_particles):
                ii_list = []
                for jj in range(ii+1, self.num_particles):
                    distances = self.positions[ii] - self.positions[jj]
                    distances = distances - (np.round(distances/self.box_size)) * self.box_size
                    distance = np.sqrt(np.sum(distances**2))
                    if distance < self.cutoff:
                        ii_list.append(jj)
                
                self.neighbours.append(ii_list.copy())

    def compute_forces(self):
        """Shifted forces using LJ potential."""

        adj = (4*self.epsilon/self.maxDist) * ((12*(self.sigma/self.maxDist)**12)-(6*(self.sigma/self.maxDist)**6))
        self.forces = np.zeros((self.num_particles, 2))  # Reset forces
        for ii in range(self.num_particles):
            for jj in self.neighbours[ii]:
                distances = self.positions[ii] - self.positions[jj]
                distances = distances - (np.round(distances/self.box_size)) * self.box_size
                distance = np.sqrt(np.sum(distances**2))
                LJforce = ((4*self.epsilon/distance) * ((12*(self.sigma/distance)**12)-(6*(self.sigma/distance)**6))) - adj
                self.forces[ii] = self.forces[ii] + ((distances/distance) * LJforce)
                self.forces[jj] = self.forces[jj] - ((distances/distance) * LJforce)

    def langevin_force(self):
        """Compute stochastic white noise and friction forces."""
        noise = np.sqrt(2 * self.kB * self.temperature * self.gamma / self.dt) * np.random.randn(self.num_particles, 2)
        return -self.gamma * self.velocities + noise

    def velocity_verlet_nve(self):
        """Velocity Verlet integration for NVE dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        self.apply_pbc()
        self.compute_forces()
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def velocity_verlet_langevin(self):
        """Velocity Verlet integration for Langevin dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        self.apply_pbc()
        self.compute_forces()
        self.forces += self.langevin_force()
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def compute_temperature(self):
        """Compute the temperature of the system from the kinetic energy."""
        kinetic_energy = 0.5 * self.mass * np.sum(self.velocities ** 2)
        return kinetic_energy / self.num_particles

    def compute_potentialenergy(self):

        potential_energy = 0
        adj = 4 * self.epsilon * ((self.sigma/self.maxDist)**12-(self.sigma/self.maxDist)**6)
        adjDer = (- 4 * self.epsilon * (12*((self.sigma/self.maxDist)**12)-6*((self.sigma/self.maxDist)**6))) / self.maxDist
        for ii in range(self.num_particles):
            for jj in self.neighbours[ii]:
                distances = self.positions[ii] - self.positions[jj]
                distances = distances - (np.round(distances/self.box_size)) * self.box_size
                distance = np.sqrt(np.sum(distances**2))
                potential_energy += (4*self.epsilon*((self.sigma/distance)**12-(self.sigma/distance)**6) - adj - ((distance - self.maxDist)*adjDer))

        return potential_energy
    
    def compute_kineticenergy(self):
        return (0.5 * self.mass * np.sum(self.velocities ** 2))

    def compute_msd(self):
        """Compute the mean squared displacement."""
        displacement = self.positions - self.initial_positions
        msd = np.mean(np.sum(displacement ** 2, axis=1))
        return msd

def part_evolution(num_particles, positions, md, points, step):
    """Animation of the particles."""

    fig = plt.figure()
    plt.xticks(np.arange(-md.box_size[0]/2, md.box_size[0]/2 + md.box_size[0]/md.division, step=md.box_size[0]/md.division))
    plt.yticks(np.arange(-md.box_size[1]/2, md.box_size[1]/2 + md.box_size[1]/md.division, step=md.box_size[1]/md.division))
    plt.grid(color = 'lightgrey', linestyle = '--', linewidth = 0.5)
    scat_dic = {}
    line_dic = {}
    for ii in range(num_particles):
        scat_dic["scat{0}".format(ii)] = plt.scatter(positions[ii, 0, 0], positions[ii, 1, 0], s=15, color='steelblue')
        #line_dic["line{0}".format(ii)] = plt.plot(positions[ii, 0, 0], positions[ii, 1, 0])[0]       
    plt.xlim([-md.box_size[0]/2, md.box_size[0]/2])
    plt.ylim([-md.box_size[1]/2, md.box_size[1]/2])
    plt.title("N=%d" %(md.num_particles))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal')

    # Only take some positions, always keep the last one
    updated_positions = positions[:, :, ::step]
    temp = np.reshape(positions[:, :, -1], (positions.shape[0], positions.shape[1], 1))
    updated_positions = np.concatenate((updated_positions, temp), axis = 2)

    axtext = fig.add_axes([0.44,0.96,0,0])
    axtext.axis("off")
    time = axtext.text(0.5,0.5, str(0), ha="left", va="top")
    
    def update(frame):
        dic = {}
        data = {}
        past = 0
        past = (0 if frame < points else (frame-points))

        time.set_text("Frame: "+str(frame)+"/"+str(updated_positions.shape[2]))

        for ii in range(num_particles):
            dic["x{0}".format(ii)] = updated_positions[ii, 0, past:frame]
            dic["y{0}".format(ii)] = updated_positions[ii, 1, past:frame]
            data["data{0}".format(ii)] = np.stack([dic["x{0}".format(ii)], dic["y{0}".format(ii)]]).T
            scat_dic["scat{0}".format(ii)].set_offsets(data["data{0}".format(ii)])
            #line_dic["line{0}".format(ii)].set_xdata(dic["x{0}".format(ii)])
            #line_dic["line{0}".format(ii)].set_ydata(dic["y{0}".format(ii)])

        return ((scat_dic["scat{0}".format(ii)]) for ii in range(num_particles))
        #return ((line_dic["line{0}".format(ii)]) for ii in range(num_particles))

    ani = animation.FuncAnimation(fig = fig, func = update, frames = updated_positions.shape[2], blit=True)
    #ani.save('images/animation.gif', writer='imagemagick', fps=30)
    plt.show()

if __name__ == '__main__':
    start = time.time()
    # Read input parameters
    directory = sys.argv[1] # Directory for input and output
    integrator = sys.argv[2] # Integrator type - options are NVE and Langevin
    num_particles = int(sys.argv[3])
    temperature = float(sys.argv[4])
    num_steps = int(float(sys.argv[5])) # Number of integration steps
    save_freq = int(num_steps/100)
    print_freq = int(num_steps/10)
    neighbour_update = 5
    
    # Create md object with input settings - more settings can be added
    md = MolecularDynamics(num_particles, temperature)

    # Create arrays for storing energy and msd
    temp = np.empty(0)
    msd = np.empty(0)
    potential = np.empty(0)
    kinetic = np.empty(0)
    total = np.zeros((md.initial_positions.shape[0], md.initial_positions.shape[1], num_steps + save_freq))
    total[:, :, 0] = md.initial_positions
    # Run integration, store and print data at given frequency
    for step in range(num_steps + save_freq):
        if step % neighbour_update == 0:
            md.compute_disk_neighbours()
            #md.compute_cell_neighbours()
        if integrator == 'nve':
            md.velocity_verlet_nve()
            total[:, :, step] = md.positions
        elif integrator == 'langevin':
            md.velocity_verlet_langevin()
        if step % save_freq == 0:
            temp = np.append(temp, md.compute_temperature())
            msd = np.append(msd, md.compute_msd())
            potential = np.append(potential, md.compute_potentialenergy())
            kinetic = np.append(kinetic, md.compute_kineticenergy())
        if step % print_freq == 0:    
            print(f"Step {step}: Kinetic = {kinetic[-1]:.4f}, Total = {potential[-1]+kinetic[-1]:.4f}")

    print("It took %fs" %(time.time()-start))
    # Plot in a gif the particles moving
    # part_evolution(num_particles, total, md, 1, 100)
    
    # Store time, temperature and energy in a single file
    time = np.arange(0, num_steps + save_freq, save_freq) * md.dt # Define time array
    np.savetxt(directory + os.sep + 'md_data.dat', np.column_stack((time, temp, potential+kinetic)))

    # Plot energy versus time
    fig, ax = plt.subplots(3, 1, figsize = (7, 7), sharex = True, dpi = 120)
    ax[0].plot(time, kinetic, color='k', linestyle='solid', marker='o', markersize='4', fillstyle='none')
    ax[0].tick_params(axis='both', labelsize=14)
    ax[0].set_ylabel("$Kinetic energy,$ $K$", fontsize=16)
    ax[1].plot(time, potential, color='k', linewidth=0.9, linestyle='solid', marker='o', markersize='6', fillstyle='none')
    ax[1].tick_params(axis='both', labelsize=14)
    ax[1].set_ylabel("$Potential Energy,$ $U$", fontsize=16)
    ax[2].plot(time, potential+kinetic, color='k', linewidth=0.9, linestyle='solid', marker='o', markersize='6', fillstyle='none')
    ax[2].tick_params(axis='both', labelsize=14)
    ax[2].set_ylabel("$Total Energy,$ $E_{tot}$", fontsize=16)
    ax[2].set_xlabel("$Simulation$ $time,$ $t$", fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.savefig("/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test/energies.png", transparent=False, format="png")

    # Plotting the potential, force and cutoff
    fig, ax = plt.subplots(2, 1, figsize = (7, 7), sharex = True, dpi = 120)
    dist = np.linspace(md.sigma*0.99, md.cutoff, 1000)
    ax[0].axhline(y=0, color="gray", linestyle="--")
    ax[0].axvline(x=md.cutoff, color="gray", linestyle="--")
    ax[0].plot(dist, 4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6), label="LJ potential")
    ax[0].plot(dist, (4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6)) - (4*md.epsilon*((md.sigma/md.cutoff)**12-(md.sigma/md.cutoff)**6)) + (dist-md.cutoff)*((4*md.epsilon/md.cutoff) * ((12*(md.sigma/md.cutoff)**12)-(6*(md.sigma/md.cutoff)**6))), label="Shifted LJ potential")
    ax[0].legend()
    ax[1].axhline(y=0, color="gray", linestyle="--")
    ax[1].axvline(x=md.cutoff, color="gray", linestyle="--")
    ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6))), label="LJ force")
    ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6)))-((4*md.epsilon/md.cutoff) * ((12*(md.sigma/md.cutoff)**12)-(6*(md.sigma/md.cutoff)**6))), label="Shifted LJ force")
    ax[1].set_ylim(top=3)
    ax[1].legend()
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.savefig("/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test/potential.png", transparent=False, format="png")
