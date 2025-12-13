import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import curve_fit
from numba import jit, njit, prange
import time
#np.random.seed(0)
#python md.py '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test' 'nve' 'WCA' 100000

def fitFunc_pow(x, a, b, c):
    return a * (x**b) + c

def fitFunc_lin(x, a, b):
    return a * x + b

def fitFunc_exp(x, a, b, c, d):
    return a * np.exp(-(x**b)/c) + d

def smallOrder():
    print("-------------------------------------------------")

def bigOrder():
    print("=================================================")

#----------------------------------------------CLASS----------------------------------------------------

class MolecularDynamics:

    def __init__(self, num_particles=100, temperature=0.1, gamma=0.1, potentialType="WCA", dt=0.0001, 
                 mass=1.0, Lx=30, Ly=30, cutoff=3, initialConf=False, figures=False, interaction=True):

        self.num_particles = num_particles
        self.potentialType = potentialType
        self.mass = mass
        self.dt = dt
        self.gamma = gamma # Friction coefficient for Langevin dynamics
        self.temperature = temperature
        self.kB = 1.0  # Boltzmann constant (arbitrary units)
        self.box_size = np.array([Lx, Ly])
        self.neighbours = [] # Initializing neighbour list
        self.sigma = 1.0  # Sigma value for potential
        if self.potentialType == "WCA":
            self.cutoff = (2**(1/6))*self.sigma # Potential cutoff
        elif self.potentialType == "LJ":
            self.cutoff = cutoff # Potential cutoff
        self.skin = 0.3 * self.sigma
        self.cellDivision = int(np.floor(self.box_size[0]/(self.cutoff + self.skin)))  # How much to devide the total space for cell neighbours
        self.epsilon = 1.0  # Epsilon value for LJ potential
        bigOrder()
        print(f"Created md object with settings:")
        print(f"Number of particles: {self.num_particles:d}\nTemperature: {self.temperature:.4f}")
        print(f"Time step: {self.dt:.4f}\nBox size: Lx {Lx:.4f} and Ly {Ly:.4f}")
        if interaction : print(f"Potential type: {self.potentialType}")
        else : print("Free particles")
        data = np.loadtxt(directory + os.sep + 'md_conf.dat')
        
        if initialConf:
            # Take the existing saved configuration
            data = np.loadtxt(directory + os.sep + 'md_conf.dat')
            self.positions = data[:, :2]
            self.velocities = data[:, 2:]
        else:
            # Initialize positions in a grid 
            nSide = int(np.ceil(np.sqrt(self.num_particles)))
            # Create grid positions
            x = (np.arange(nSide) + 0.5) * (self.box_size[0]/nSide) - (self.box_size[0]/nSide)/2
            y = (np.arange(nSide) + 0.5) * (self.box_size[1]/nSide) - (self.box_size[1]/nSide)/2
            xv, yv = np.meshgrid(x, y)
            positions = np.vstack([xv.ravel(), yv.ravel()]).T
            positions = positions - (self.box_size[0]/2)*np.ones_like(positions)
            positions = (positions + self.box_size / 2) % self.box_size - self.box_size / 2
            self.positions = positions[:self.num_particles] # Only return the first num_particles if grid has extra points

            # Initialize velocities
            self.velocities = np.random.normal(0, np.sqrt((self.kB*self.temperature)/self.mass), (self.num_particles, 2)) # Maxwell-Boltmann
            #self.velocities = self.velocities - (np.sum(self.velocities, axis=0)/self.num_particles) # Remove center of mass
            """# Re-sample outliers
            vMax = self.cutoff/3
            outliers = np.sqrt(np.sum(self.velocities**2, axis=1)) > vMax
            while outliers.any():
                self.velocities[outliers] = np.random.normal(0, 1, (np.sum(outliers), 2))
                self.velocities = self.velocities - (np.sum(self.velocities, axis=0)/self.num_particles) # Remove center of mass
                outliers = np.sqrt(np.sum(self.velocities**2, axis=1)) > vMax"""
            self.velocities = self.velocities * np.sqrt((self.num_particles * self.temperature)/(0.5 * self.mass * np.sum(self.velocities ** 2))) # Scale to have initial temperature
        
        self.initial_positions = np.copy(self.positions) # Store initial positions to compute the MSD
        self.lastsaved_positions = np.copy(self.positions) # Store the last needed configuration for the neighbours update
        self.forces = np.zeros((num_particles, 2))
        self.forcesContainer = []
        self.allforcesContainer = []
        self.distancesContainer = []
        self.potentialEnergy = 0
        self.figures = figures
        self.unwrapped_positions = np.copy(self.positions)
        self.all_positions = []
        self.all_unwrapped_positions = []
        self.positions_save_freq = 1000

        # Set size for interacting particles - no force is implemented yet
        self.interaction = interaction
        """if self.interaction:
            # Set particle sizes based on the shortest interparticle distance.
            min_distance = np.inf
            for i in range(self.num_particles):
                for j in range(i + 1, self.num_particles):
                    r_ij = np.linalg.norm(self.positions[i] - self.positions[j])
                    if r_ij < min_distance:
                        min_distance = r_ij
            # Set all particles' radii to the minimum interparticle distance
            self.radii = np.full(self.num_particles, min_distance)
            print(f'Average particle diameter: {np.mean(self.radii).df}')"""

    def apply_pbc(self):
        """Apply periodic boundary conditions to keep particles inside the simulation box."""
        self.positions = (self.positions + self.box_size / 2) % self.box_size - self.box_size / 2

    def compute_cell_neighbours(self):
        """"Computing nearest neighbours based on cell subdivision."""

        self.neighbours = [] # Reset neighbours
        head = -np.ones((self.cellDivision, self.cellDivision), dtype=int)
        cell = np.zeros((self.num_particles, 2), dtype=int)
        list = np.zeros(self.num_particles, dtype=int)
        near_cells = [(1, 1), (1, 0), (1, -1), (0, 1)]
        other_cell = [0, 0]

        for ii in range(self.num_particles):
            cell[ii, :] = np.floor(((self.positions[ii, :] + 5)*self.cellDivision)/(self.box_size[0])) # which cell ii belongs to
            list[ii] = head[cell[ii, 0], cell[ii, 1]] # point ii to previous head of the cell, 0 if it is first in cell
            head[cell[ii, 0], cell[ii, 1]] = ii # ii is now the new head of the cell

        for ii in range(self.num_particles):
            ii_list = []
            current_cell = (cell[ii, 0], cell[ii, 1])
            current_head = head[cell[ii, 0], cell[ii, 1]]

            # Neighbours in current cell
            while current_head != -1: # if there are other neighbours and the cell wasn't checked...
                if current_head > ii: # ...and it's not yourself (and no double count)
                    ii_list.append(current_head) # append neighbour
                current_head = list[current_head] # next head in line

            # Neighbours in near cells
            for value in near_cells:
                other_cell[0] = current_cell[0] + value[0]
                other_cell[1] = current_cell[1] + value[1]
                other_head = head[other_cell[0]%self.cellDivision, other_cell[1]%self.cellDivision]
                while other_head != -1: # if there are other neighbours...
                    ii_list.append(other_head) # append neighbour
                    other_head = list[other_head] # next head in line
                    
            self.neighbours.append(ii_list.copy())
        
        self.lastsaved_positions = self.positions

    def compute_disk_neighbours(self):
        """"Computing nearest neighbours based on simple disk distance."""

        self.neighbours = [] # Reset neighbours

        for ii in range(self.num_particles):
                ii_list = []
                for jj in range(ii+1, self.num_particles):
                    distances = self.positions[ii] - self.positions[jj]
                    distances = distances - (np.round(distances/self.box_size)) * self.box_size
                    distance = np.sqrt(np.sum(distances**2))
                    if distance <= self.cutoff + self.skin:
                        ii_list.append(jj)
                
                self.neighbours.append(ii_list.copy())

        self.lastsaved_positions = self.positions

    def compute_LJ_forces(self):
        """Shifted forces using LJ potential."""
        
        forceShift = (4*self.epsilon/self.cutoff) * ((12*(self.sigma/self.cutoff)**12)-(6*(self.sigma/self.cutoff)**6))
        potShift = 4 * self.epsilon * ((self.sigma/self.cutoff)**12-(self.sigma/self.cutoff)**6)
        potDerShift = (- 4 * self.epsilon * (12*((self.sigma/self.cutoff)**12)-6*((self.sigma/self.cutoff)**6))) / self.cutoff
        self.forces = np.zeros((self.num_particles, 2))  # Reset forces
        potential_energy = 0  # Reset potential
        for ii in range(self.num_particles):
            for jj in self.neighbours[ii]:
            #for jj in range(ii + 1, self.num_particles):
                distances = self.positions[ii] - self.positions[jj]
                distances = distances - (np.round(distances/self.box_size)) * self.box_size
                distance = np.sqrt(np.sum(distances**2))
                if distance < self.cutoff:
                    potential_energy += (4*self.epsilon*((self.sigma/distance)**12-(self.sigma/distance)**6) - potShift - ((distance - self.cutoff)*potDerShift))
                    LJforce = ((4*self.epsilon/distance) * ((12*(self.sigma/distance)**12)-(6*(self.sigma/distance)**6))) - forceShift
                    self.forces[ii] = self.forces[ii] + ((distances/distance) * LJforce)
                    self.forces[jj] = self.forces[jj] - ((distances/distance) * LJforce)
                    #self.forcesContainer.append(LJforce)
                    #self.distancesContainer.append(distance)

        self.potentialEnergy = potential_energy

    def compute_WCA_forces(self):
        """Forces using WCA potential."""
        
        self.forces = np.zeros((self.num_particles, 2))  # Reset forces
        potential_energy = np.zeros(self.num_particles)  # Reset potential
        for ii in range(self.num_particles):
            for jj in self.neighbours[ii]:
            #for jj in range(ii + 1, self.num_particles):
                distances = self.positions[ii] - self.positions[jj]
                distances -= np.round(distances/self.box_size) * self.box_size
                distance = np.linalg.norm(distances)
                WCAforce = 0
                ratio6 = (self.sigma / distance)**6
                ratio12 = ratio6 * ratio6
                if distance < self.cutoff:
                    # 0.5 for distributing the energy in the two particles
                    potential_energy[ii] += 0.5 * self.epsilon * (4 * (ratio12 - ratio6) + 1)
                    potential_energy[jj] += 0.5 * self.epsilon * (4 * (ratio12 - ratio6) + 1)
                    #potential_energy[ii] += 0.5 * (self.epsilon * 4 * ((self.sigma/distance)**12 - (self.sigma/distance)**6) + self.epsilon)
                    #potential_energy[jj] += 0.5 * (self.epsilon * 4 * ((self.sigma/distance)**12 - (self.sigma/distance)**6) + self.epsilon)
                    #WCAforce = (4 * self.epsilon / distance) * (12 * (self.sigma/distance)**12 - 6 * (self.sigma/distance)**6)
                    WCAforce = 24 * self.epsilon * (2 * ratio12 - ratio6) / distance 
                    self.forces[ii] += WCAforce * distances / distance
                    self.forces[jj] -= WCAforce * distances / distance
                #if (ii==2) and (jj==5):
                #    self.allforcesContainer.append(WCAforce)
                #if abs(distance - self.cutoff) < 0.01:
                #    self.forcesContainer.append(WCAforce)
                #    self.distancesContainer.append(distance)

        self.potentialEnergy = np.sum(potential_energy)

    def langevin_force(self):
        """Compute stochastic white noise and friction forces."""
        noise = np.sqrt(2 * self.kB * self.temperature * self.gamma / self.dt) * np.random.randn(self.num_particles, 2)
        return -self.gamma * self.velocities + noise

    def velocity_verlet_nve(self):
        """Velocity Verlet integration for NVE dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        self.unwrapped_positions += self.velocities * self.dt
        self.apply_pbc()
        if self.interaction:
            if self.potentialType == "WCA":
                self.compute_WCA_forces()
            elif self.potentialType == "LJ":
                self.compute_LJ_forces()
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def velocity_verlet_langevin(self):
        """Velocity Verlet integration for Langevin dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        self.unwrapped_positions += self.velocities * self.dt
        self.apply_pbc()
        if self.interaction:
            if self.potentialType == "WCA":
                self.compute_WCA_forces()
            elif self.potentialType == "LJ":
                self.compute_LJ_forces()
        else : self.forces = np.zeros((num_particles, 2))
        self.forces += self.langevin_force()
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def compute_temperature(self):
        """Compute the temperature of the system from the kinetic energy."""
        kinetic_energy = 0.5 * self.mass * np.sum(self.velocities ** 2)
        return kinetic_energy / self.num_particles
    
    def compute_kineticenergy(self):
        return (0.5 * self.mass * np.sum(self.velocities ** 2))

    def compute_msd(self):
        """Compute the Mean Squared Displacement."""
        displacement = self.unwrapped_positions - self.initial_positions
        msd = np.mean(np.sum(displacement ** 2, axis=1))
        return msd

    def compute_isf(self, inf_lim, sup_lim):
        "Compute the Intermediate Scattering Function for different k modulus with mean over different t0."

        # print("Computing ISF with t every", self.positions_save_freq, "steps.")
        angles = np.arange(0, 2*np.pi, np.pi/4)
        isf_self = np.zeros((np.shape(self.kmods)[0]), dtype=complex)
        isf_int = np.zeros((np.shape(self.kmods)[0]), dtype=complex)
        #mask = ~np.eye(self.num_particles, dtype=bool)
        isf_total_self = np.zeros((sup_lim, (np.shape(self.kmods)[0])), dtype=complex)
        isf_total_int = np.zeros((sup_lim, (np.shape(self.kmods)[0])), dtype=complex)

        for t0 in range(inf_lim, sup_lim):
            temporary_ip = self.all_unwrapped_positions[t0]
            for t in range(t0, sup_lim):
                isf_self = np.zeros_like(isf_self)
                isf_int = np.zeros_like(isf_int)
                for ii, kMod in enumerate(self.kmods):
                    for angle in angles:
                        kVec = np.array((kMod * np.cos(angle), kMod * np.sin(angle)))
                        isf_self[ii] += np.sum(np.exp(1j * np.matmul((self.all_unwrapped_positions[t] - temporary_ip), kVec)))/(self.num_particles * np.shape(angles)[0])
                        isf_int[ii] += np.sum(np.exp(1j * np.matmul((self.all_unwrapped_positions[t][ :, None, :] - temporary_ip[None, :, :]), kVec)))/(self.num_particles*(self.num_particles-1) * np.shape(angles)[0])
                isf_total_self[t-t0, :] += (isf_self) / ((sup_lim - inf_lim) - (t-t0))
                isf_total_int[t-t0, :] += (isf_int) / ((sup_lim - inf_lim) - (t-t0))

        if (np.any(np.abs(np.imag(isf_total_self)) > 1e-4) or np.any(np.abs(np.imag(isf_total_int))  > 1e-4)):
            print("The imaginary part is absolutely too big.")

        isf_total_self = np.real(isf_total_self)
        isf_total_int = np.real(isf_total_int)

        return isf_total_self, isf_total_int

#-----------------------------------------OTHER-FUNCTIONS-----------------------------------------------

def part_evolution(num_particles, positions, md, points, step):
    """Animation of the particles."""

    fig = plt.figure()
    ax = fig.add_subplot(111)
    subdivision = 10
    plt.xticks(np.arange(-md.box_size[0]/2, md.box_size[0]/2 + md.box_size[0]/subdivision, step=md.box_size[0]/subdivision))
    plt.yticks(np.arange(-md.box_size[1]/2, md.box_size[1]/2 + md.box_size[1]/subdivision, step=md.box_size[1]/subdivision))
    #plt.grid(color = 'lightgrey', linestyle = '--', linewidth = 0.5)
    scat_dic = {}
    line_dic = {}

    # Represent the correct size of the diameter
    radius = md.sigma/2
    trans = ax.transData.transform
    inv = fig.dpi_scale_trans.inverted().transform  
    x0, y0 = trans((0,0))
    x1, y1 = trans((radius, 0))
    radius_pixels = x1 - x0

    for ii in range(num_particles):
        scat_dic["scat{0}".format(ii)] = plt.scatter(positions[ii, 0, 0], positions[ii, 1, 0], s=radius_pixels**2, 
                                                     facecolor='goldenrod', edgecolor="black", linewidth=0.7)
        #line_dic["line{0}".format(ii)] = plt.plot(positions[ii, 0, 0], positions[ii, 1, 0])[0]       
    plt.xlim([-md.box_size[0]/2, md.box_size[0]/2])
    plt.ylim([-md.box_size[1]/2, md.box_size[1]/2])
    plt.title(r"N=%d, T=%.1f, cutoff=%.1f$\sigma$" %(md.num_particles, md.temperature, md.cutoff/md.sigma))
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
    #ani.save('test/animation.gif', writer='imagemagick', fps=30)
    plt.show()

def reduce_vectors(vector1, vector2, eps):
    # v2 is the positions
    
    vector1 = np.asarray(vector1)
    vector2 = np.asarray(vector2)
    idx = np.argsort(vector2)
    v1_sorted = vector1[idx]
    v2_sorted = vector2[idx]
    groups = []
    current_group = [0]

    for ii in range(1, len(v2_sorted)):
        if abs(v2_sorted[ii] - v2_sorted[current_group[-1]]) <= eps:
            current_group.append(ii)
        else:
            groups.append(current_group)
            current_group = [ii]
    groups.append(current_group)

    reduced_v1 = np.array([v1_sorted[np.array(gg)].mean() for gg in groups])
    reduced_v2 = np.array([v2_sorted[np.array(gg)].mean() for gg in groups])

    return reduced_v1, reduced_v2

#-----------------------------------------------MAIN----------------------------------------------------

if __name__ == '__main__':

    start = time.time()
    # Read input parameters
    directory = sys.argv[1] # Directory for input and output
    integrator = sys.argv[2] # Integrator type - options are NVE and Langevin
    potentialType = sys.argv[3] # Options are LJ and WCA
    num_steps = int(float(sys.argv[4])) # Number of integration steps
    save_freq = int(num_steps/100)
    print_freq = int(num_steps/10)

    # Control what graphs to show 
    compute_energy = False
    compute_msd = True
    compute_isf = True
    compute_gif = False
    store_data = False
    interaction = False

    # Control the parameters (to have a single run, just put a single value for each)
    particles = np.array([100], dtype=int)
    temperatures = np.array([10], dtype=float)
    frictions = np.array([1.0], dtype=float)
    msd_total = np.zeros((np.shape(particles)[0], np.shape(temperatures)[0], np.shape(frictions)[0], int(num_steps/save_freq) + 1))
    kmods = np.array([2*np.pi]) 
    isf_total = np.zeros((np.shape(particles)[0], np.shape(temperatures)[0], np.shape(frictions)[0], 
                          int(num_steps/save_freq) + 1, np.shape(kmods)[0], 2))

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):

                # Create md object with input settings - more settings can be added
                md = MolecularDynamics(num_particles, temperature, gamma, potentialType, interaction = interaction)
                md.kmods = kmods/md.box_size[0]
                md.positions_save_freq = save_freq

                # Create arrays for storing energy and msd
                temp = np.empty(0)
                msd = np.empty(0)
                potential = np.empty(0)
                kinetic = np.empty(0)
                md.compute_disk_neighbours()
                #md.compute_cell_neighbours()

                # Run integration, store and print data at given frequency
                smallOrder()
                for step in range(num_steps + save_freq):
                    distances = md.positions - md.lastsaved_positions
                    distances -= np.round(distances/md.box_size) * md.box_size
                    distance = np.linalg.norm(distances, axis=1)
                    if np.any(distance >= md.skin/2): # Update the neighbour list only when necessary
                        md.compute_disk_neighbours() 
                        #md.compute_cell_neighbours()
                    if integrator == 'nve':
                        md.velocity_verlet_nve()
                    elif integrator == 'langevin':
                        md.velocity_verlet_langevin()
                    if step % save_freq == 0:
                        temp = np.append(temp, md.compute_temperature())
                        if compute_msd : msd = np.append(msd, md.compute_msd())
                        potential = np.append(potential, md.potentialEnergy)
                        kinetic = np.append(kinetic, md.compute_kineticenergy())
                    if step % md.positions_save_freq == 0:
                        md.all_positions.append(md.positions.copy())
                        md.all_unwrapped_positions.append(md.unwrapped_positions.copy())
                    if step % print_freq == 0:
                        print(f"Step {step}, T: {temp[-1]:.4f}, E: {potential[-1]+kinetic[-1]:.7f}")
                
                if compute_msd: msd_total[ii, jj, kk, :] = msd

                if compute_isf: 
                    isf_total[ii, jj, kk, :, :, 0], isf_total[ii, jj, kk, :, :, 1] = md.compute_isf(inf_lim = 0, sup_lim = np.shape(md.all_unwrapped_positions)[0])

                # Plot in a gif the particles moving
                if compute_gif: 
                    md.all_positions = np.array(md.all_positions)
                    md.all_positions = np.stack(md.all_positions, axis=-1)  
                    part_evolution(num_particles, md.all_positions, md, 1, 1)

    smallOrder()
    print("It took %fs" %(time.time()-start))
    bigOrder()
    
    time = np.arange(0, num_steps + save_freq, save_freq) * md.dt # Define time array
    if store_data:
        # Store time, temperature and energy in a single file
        np.savetxt(directory + os.sep + 'md_data.dat', np.column_stack((time, temp, potential, kinetic)))
        # Store also the last positions and velocities of the particles
        np.savetxt(directory + os.sep + 'md_conf.dat', np.column_stack((md.positions[:, 0], md.positions[:, 1], md.velocities[:, 0], md.velocities[:, 1])))
     
    # Plot energy versus time (valid only for last iteration for now)
    if compute_energy:
        plt.clf()
        plt.title(r"Energies for: N=%d ($\rho$=%.1f), T=%.1f" %(md.num_particles, md.num_particles/(md.box_size[0]**2), md.temperature), fontsize=16)
        plt.plot(time, kinetic, color='seagreen', linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Kinetic energy $K$")
        plt.tick_params(axis='both', labelsize=14)
        plt.plot(time, potential, color='steelblue', linewidth=0.9, linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Potential energy $U$")
        plt.tick_params(axis='both', labelsize=14)
        plt.plot(time, potential+kinetic, color='orchid', linewidth=0.9, linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Total energy $E_{tot}$")
        plt.tick_params(axis='both', labelsize=14)
        plt.ylabel("Energies", fontsize=14)
        plt.xlabel(r"Simulation time, $t$", fontsize=14)
        plt.tight_layout()
        plt.legend()
        plt.savefig(directory + "/energies.png", transparent=False, format="png")

    # Plot mean squared displacement
    if compute_msd:
        plt.clf()
        eq = md.gamma*md.dt*save_freq
        fit = True
        if integrator == 'nve': plt.title(r"MSD with interaction:%s and Langevin:False" %(interaction), fontsize=16)
        else : plt.title(r"MSD with interaction:%s and Langevin:True" %(interaction), fontsize=16)
        cmap = plt.get_cmap("viridis")  
        n_colors = int(np.shape(particles)[0]) * int(np.shape(temperatures)[0]) * int(np.shape(frictions)[0]) 
        colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
        counter = 0
        for ii, num_particles in enumerate(particles):
            for jj, temperature in enumerate(temperatures):
                for kk, gamma in enumerate(frictions):
                    if integrator == 'nve': 
                        if (interaction == False) and fit:
                            popt, pcov = curve_fit(fitFunc_pow, time, msd_total[ii, jj, kk, :]) 
                            plt.plot(time, fitFunc_pow(time, popt[0], popt[1], popt[2]), color=colors[counter], linestyle='solid', linewidth=1) 
                            plt.plot(time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$\rho$=%.2f, T=%.1f, $\propto t^{%.1f}$" %(num_particles/(md.box_size[0]**2), temperature, popt[1]))
                        else:
                            plt.plot(time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$\rho$=%.2f, T=%.1f" %(num_particles/(md.box_size[0]**2), temperature))
                    elif integrator == 'langevin': 
                        if (interaction == False) and fit:
                            # Ballistic regime
                            popt1, pcov1 = curve_fit(fitFunc_pow, time[:int(1/eq)], msd_total[ii, jj, kk, :int(1/eq)]) 
                            plt.plot(time[:int(1/eq)], fitFunc_pow(time[:int(1/eq)], popt1[0], popt1[1], popt1[2]), color=colors[counter], linestyle='solid', linewidth=1) 
                            # Diffusive regime
                            popt2, pcov2 = curve_fit(fitFunc_lin, time[int(6/eq):], msd_total[ii, jj, kk, int(6/eq):]) 
                            plt.plot(time[int(6/eq):], fitFunc_lin(time[int(6/eq):], popt2[0], popt2[1]), color=colors[counter], linestyle='solid', linewidth=1) 
                            plt.plot(time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f, $\propto t^{%.1f}\rightarrow\propto t$" %(num_particles/(md.box_size[0]**2), temperature, gamma, popt1[1]))
                        else:
                            plt.plot(time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f" %(num_particles/(md.box_size[0]**2), temperature, gamma))
                    counter += 1
        plt.ylabel(r"MSD, $\langle |r(t)-r_0|^2 \rangle$", fontsize=14)
        plt.xlabel(r"Simulation time, $t$", fontsize=14)
        plt.tight_layout()
        plt.ylim(bottom=0.1, top=1.5*np.max(msd_total))
        plt.xscale("log")
        plt.yscale("log")
        plt.legend()
        plt.savefig(directory + "/msd.png", transparent=False, format="png")

    # Plot intermediate scattering function
    if compute_isf:
        plt.clf()
        eq = md.gamma*md.dt*save_freq
        fit = True
        if integrator == 'nve': plt.title(r"ISF with interaction:%s and Langevin:False" %(interaction), fontsize=16)
        else : plt.title(r"ISF with interaction:%s and Langevin:True" %(interaction), fontsize=16)
        cmap = plt.get_cmap("viridis")  
        n_colors = int(np.shape(md.kmods)[0]) 
        colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
        taus = np.zeros((np.shape(temperatures)[0], np.shape(md.kmods)[0]))
        transparent = int(np.shape(particles)[0]) * int(np.shape(temperatures)[0]) * int(np.shape(frictions)[0]) 
        transparency = 1/transparent
        for ii, num_particles in enumerate(particles):
            for jj, temperature in enumerate(temperatures):
                for kk, gamma in enumerate(frictions):
                    for mm, kMod in enumerate(md.kmods):
                        #if np.shape(temperatures)[0] > 1: taus[ii, mm] = popt[1]
                        if integrator == 'nve': 
                            if (interaction == False) and fit:
                                popt, pcov = curve_fit(fitFunc_exp, time[0:np.shape(isf_total)[3]], isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1], 
                                                maxfev=100000, p0=[1, 2, 1, 0])
                                plt.plot(time[0:np.shape(isf_total)[3]], fitFunc_exp(time[0:np.shape(isf_total)[3]], popt[0], popt[1], popt[2], popt[3]), 
                                        color=colors[mm], linewidth=1, linestyle='solid', alpha = transparency)
                                plt.plot(time[0:np.shape(isf_total)[3]], isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1], 
                                         color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$\rho$=%.2f, T=%.1f, $|k|=%.1f$, $\propto e^{-(t-t_0)^{%.1f}/%.1f}$" %(num_particles/(md.box_size[0]**2), temperature, kMod, popt[1], popt[2]))
                            else:
                                plt.plot(time[0:np.shape(isf_total)[3]], isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1], 
                                         color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$\rho$=%.2f, T=%.1f, $|k|=%.1f$" %(num_particles/(md.box_size[0]**2), temperature, kMod))
                        elif integrator == 'langevin': 
                            if (interaction == False) and fit:
                                popt1, pcov1 = curve_fit(fitFunc_exp, time[:int(1/eq)], isf_total[ii, jj, kk, :int(1/eq), mm, 0] + isf_total[ii, jj, kk, :int(1/eq), mm, 1], 
                                               maxfev=100000, p0=[1, 2, 1, 0]) 
                                plt.plot(time[:int(1/eq)], fitFunc_exp(time[:int(1/eq)], popt1[0], popt1[1], popt1[2], popt1[3]), 
                                        color=colors[mm], linewidth=1, linestyle='solid', alpha = transparency)
                                # Diffusive regime
                                #popt2, pcov2 = curve_fit(fitFunc_exp, time[int(8/eq):], isf_total[ii, jj, kk, int(8/eq):, mm, 0] + isf_total[ii, jj, kk, int(8/eq):, mm, 1], 
                                #                maxfev=100000, p0=[1, 1, 1, 0]) 
                                #plt.plot(time[int(8/eq):], fitFunc_exp(time[int(8/eq):], popt2[0], popt2[1], popt2[2], popt2[3]), 
                                #        color=colors[mm], linestyle='solid', alpha = transparency)
                                plt.plot(time[0:np.shape(isf_total)[3]], isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1],
                                        color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f, $|k|=%.1f$, $\propto e^{-(t-t_0)^{%.1f}/%.1f}$" %(num_particles/(md.box_size[0]**2), temperature, gamma, kMod, popt1[1], popt1[2]))
                            else:
                                plt.plot(time[0:np.shape(isf_total)[3]], isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1],
                                         color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f, $|k|=%.1f$" %(num_particles/(md.box_size[0]**2), temperature, gamma, kMod))
                    transparency += 1/transparent
        plt.ylabel(r"ISF", fontsize=14)
        plt.xlabel(r"Simulation time, $t-t_0$", fontsize=14)
        #plt.ylim(bottom=0)
        #plt.xscale("log")
        plt.axhline(y=0, color="gray", linestyle="--")
        plt.tight_layout()
        plt.legend()
        plt.savefig(directory + "/isf.png", transparent=False, format="png")

        """if np.shape(temperatures)[0] > 1:
            plt.clf()
            if integrator == 'nve': plt.title(r"$\tau$ with interaction:%s, increasing $T$" %(interaction), fontsize=16)
            else : plt.title(r"\tau$ with interaction:%s, increasing $T and Langevin" %(interaction), fontsize=16)
            for mm, kMod in enumerate(md.kmods):
                plt.plot(temperatures, taus[:, mm], color=colors[mm], linestyle='solid', marker='o', markersize='4', fillstyle='none', label=r"$|k|=%.1f$" %(kMod))
            plt.ylabel(r"$\tau$", fontsize=14)
            plt.xlabel(r"Temperature $T$", fontsize=14)
            plt.tight_layout()
            plt.legend()
            plt.savefig(directory + "/tau.png", transparent=False, format="png")"""

    # Other studies
    if md.figures:
        # Plotting the potential, force and cutoff
        fig, ax = plt.subplots(2, 1, figsize = (7, 7), sharex = True, dpi = 120)
        dist = np.linspace(md.sigma*0.9, md.cutoff, 1000)
        added = np.linspace(md.cutoff, md.cutoff*1.1, 1000)
        #forces, distances = reduce_vectors(md.forcesContainer, md.distancesContainer, 1e-06)
        ax[0].axhline(y=0, color="gray", linestyle="--")
        ax[0].axvline(x=md.cutoff, color="gray", linestyle="--")
        ax[1].axhline(y=0, color="gray", linestyle="--")
        ax[1].axvline(x=md.cutoff, color="gray", linestyle="--")
        if md.potentialType == "WCA":
            ax[0].plot(dist, 4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6) + md.epsilon, label="WCA potential", color="orchid")
            ax[0].plot(added, 0*added, color="orchid")
            ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6))), label="WCA force", color="orchid")
            ax[1].plot(added, 0*added, color="orchid")
        elif md.potentialType == "LJ":
            ax[0].plot(dist, 4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6), label="LJ potential")
            ax[0].plot(dist, (4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6)) - (4*md.epsilon*((md.sigma/md.cutoff)**12-(md.sigma/md.cutoff)**6)) + (dist-md.cutoff)*((4*md.epsilon/md.cutoff) * ((12*(md.sigma/md.cutoff)**12)-(6*(md.sigma/md.cutoff)**6))), label="Shifted LJ potential")
            ax[0].plot(added, 0*added, color="orange")
            ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6))), label="LJ force")
            ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6)))-((4*md.epsilon/md.cutoff) * ((12*(md.sigma/md.cutoff)**12)-(6*(md.sigma/md.cutoff)**6))), label="Shifted LJ force")
            ax[1].plot(added, 0*added, color="orange")
        ax[0].legend()
        ax[1].set_ylim(top=100)
        ax[1].scatter(md.distancesContainer, md.forcesContainer, color="lime", s=1, label="All forces")
        ax[1].legend()
        plt.tight_layout()
        plt.subplots_adjust(hspace=0)
        plt.savefig(directory + "/potential.png", transparent=False, format="png")

        # Other checks
        plt.clf()
        plt.axhline(y=0, color="gray", linestyle="--")
        plt.axvline(x=md.cutoff, color="gray", linestyle="--")
        plt.scatter(md.distancesContainer, md.forcesContainer, color="orchid", s=4, label="Zoomed forces (all pairs)")
        plt.legend()
        plt.savefig(directory + "/continuous.png", transparent=False, format="png")

        # Other checks
        plt.clf()
        plt.axhline(y=0, color="gray", linestyle="--")
        plt.plot(md.allforcesContainer, color="orchid", label="Zoomed forces over time (sampled pair)")
        plt.ylim(top=0.01, bottom=-0.001)
        plt.legend()
        plt.savefig(directory + "/continuous2.png", transparent=False, format="png")
