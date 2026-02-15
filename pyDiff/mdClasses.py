"""
mdClasses 
----------
Here are contained the used classes to simulate molecular dynamics, mainly,
the class MolecularDynamics.
"""
import os
import numpy as np
import numba as nb
from mdFunctions import *
home = '/home/auroisflying/thesis/simSoft/pyDiff/test'

@nb.njit(parallel=True, fastmath=True)
def compute_WCA_forces_numba(positions, neighbours, neighbour_counts, box_size, sigma, epsilon, cutoff):
    num_particles = positions.shape[0]
    forces = np.zeros((num_particles, 2))
    energy = np.zeros(num_particles)  # Reset potential

    for ii in nb.prange(num_particles):
        for k in range(neighbour_counts[ii]):
            jj = neighbours[ii, k]
            #if ii == 0: print("particle ii:", ii, "neighbor:", jj)
            #if jj == 0: print("neighbor jj:", jj, "particle:", ii)
            distances = positions[ii] - positions[jj]
            distances -= np.round(distances/box_size) * box_size
            distance = np.linalg.norm(distances)
            if distance < cutoff:
                ratio6 = (sigma / distance)**6
                ratio12 = ratio6 * ratio6
                energy[ii] += 0.5 * 0.5 * epsilon * (4 * (ratio12 - ratio6) + 1)
                WCAforce = 24 * epsilon * (2 * ratio12 - ratio6) / distance 
                forces[ii] += 0.5 * WCAforce * distances / distance

    return forces, energy

@nb.njit(parallel=True, fastmath=True)
def compute_disk_neighbours_numba(num_particles, positions, box_size, cutoff, skin, neighbours, max_neighbors):
        """"Computes nearest neighbours based on disk distance."""
        
        neighbour_counts = np.zeros(num_particles, dtype=np.int32)

        for ii in nb.prange(num_particles):
            count = 0
            for jj in range(num_particles):
                if jj != ii:
                    distances = positions[ii] - positions[jj]
                    distances -= np.round(distances / box_size) * box_size
                    if np.sum(distances**2) <= (cutoff + skin)**2:
                        if count < max_neighbors:
                            neighbours[ii, count] = jj
                            count += 1
            neighbour_counts[ii] = count

        neighborCheckPositions = positions
        return neighborCheckPositions, neighbour_counts, neighbours

class MolecularDynamics:

    def __init__(self, num_particles: int = 100, temperature: float = 1.0, gamma: float = 1.0, potentialType: str = "WCA", integrator: str = 'nve',
                 dt: float = 0.0001, Lx: int = 10, Ly: int = 10, initialConf: bool = False, interaction: bool = False, mixture = False, steps = 1e03):
        """
        This class implements a series of methods to simulate molecular dynamics with different potential types and
        different algorithms in 2D.

        Parameters
        ----------
        num_particles : int
            Number of particles of the system, default is 100.
        temperature : float
            Temperature of the system, default is 1.0.
        gamma : float
            Friction of the system, default is 1.0.
        potentialType : str
            What potential to use. Options are "WCA" and "LJ", default is "WCA".
        integrator : str
            The integrator used for the object. Options are 'nve', 'langevin and 'em', default is 'nve'.
        dt : float
            Timestep of the simulation, default id 0.0001.
        Lx, Ly : int
            Dimensions of the 2D box.
        initialConf : bool
            If True, loads an existing configuration. Default is False.
        interaction : bool
            If True, there is interaction. Default is False.
        """

        # Variables set to specific values
        self.mass = 1.0
        self.kB = 1.0 
        self.sigma = 1.0 
        self.epsilon = 1.0

        # Initialization from input
        self.num_particles = num_particles
        self.dt = dt 
        self.gamma = gamma 
        self.temperature = temperature
        self.box_size = np.array([Lx, Ly])
        self.interaction = interaction
        self.mixture = mixture
        self.potentialType = potentialType
        self.integrator = integrator
        self.steps = steps

        # Potential
        if self.potentialType == "WCA" or self.potentialType == "WCAnumba":
            self.cutoff = (2**(1/6))*self.sigma 
        elif self.potentialType == "LJ":
            self.cutoff = 3.0 
        else:
            self.cutoff = 2.0
        self.forces = np.zeros((num_particles, 2))
        self.energy = np.zeros(num_particles) # Per-perticle potential energy array
        self.potentialEnergy = 0
        # Neighbours list
        self.max_neighbors = 64  # Choose safely (depends on density)
        self.neighbours = np.full((self.num_particles, self.max_neighbors), -1, dtype=np.int32)
        self.neighbour_counts = np.zeros(self.num_particles, dtype=np.int32)
        self.skin = 0.3 * self.sigma 
        self.cellDivision = int(np.floor(max(self.box_size[0], self.box_size[1])/(self.cutoff + self.skin))) 
        # Positions
        if initialConf:
            # Take the existing saved configuration from the first iteration
            if self.interaction : optionsDirectory = f"{self.integrator}WCA_N{self.num_particles:d}_phi{self.num_particles * (np.pi * (0.5)**2)/(self.box_size[0]*self.box_size[1]):.1f}_T{self.temperature:.1f}_g{self.gamma:.2f}"
            else : optionsDirectory = f"{integrator}FREE_N{self.num_particles:d}_T{self.temperature:.1f}_g{self.gamma:.2f}"
            savePath = os.path.join(home, optionsDirectory)
            data = np.load(savePath + os.sep + 'initialConfiguration.npz')
            self.positions = data["positions"]
            #self.velocities = data[:, 2:]
        else:
            # Initialize positions in a grid 
            nx = int(np.ceil(np.sqrt(self.num_particles * self.box_size[0] / self.box_size[1])))
            ny = int(np.ceil(self.num_particles / nx))
            # Create grid positions
            x = (np.arange(nx) + 0.5) * (self.box_size[0]/nx) - (self.box_size[0]/nx)/2
            y = (np.arange(ny) + 0.5) * (self.box_size[1]/ny) - (self.box_size[1]/ny)/2
            xv, yv = np.meshgrid(x, y)
            positions = np.vstack([xv.ravel(), yv.ravel()]).T
            positions[:, 0] = positions[:, 0] - (self.box_size[0]/2)*np.ones_like(positions[:, 0])
            positions[:, 1] = positions[:, 1] - (self.box_size[1]/2)*np.ones_like(positions[:, 1])
            positions = (positions + self.box_size / 2) % self.box_size - self.box_size / 2
            self.positions = positions[:self.num_particles] # Only return the first num_particles if grid has extra points

        self.neighborCheckPositions = np.copy(self.positions)
        self.unwrappedPositions = np.copy(self.positions)
        self.allPositions = []
        self.allUnwrappedPositions = []
        self.allVelocities = []
        # Velocities (Maxwell-Boltzmann without center of mass and scaled to correct initial temperature)
        self.velocities = np.random.normal(0, np.sqrt((self.kB*self.temperature)/self.mass), (self.num_particles, 2)) 
        self.velocities = self.velocities - (np.sum(self.velocities, axis=0)/self.num_particles) 
        self.velocities = self.velocities * np.sqrt((self.num_particles * self.temperature)/(0.5 * self.mass * np.sum(self.velocities ** 2))) 
        #print("Center of mass velocity: ", np.sum(self.velocities, axis=0)/self.num_particles)
        # Activity variables
        self.tau = 10
        #self.activityForce = np.sqrt(2 * self.kB * self.temperature * self.gamma / self.tau)
        self.thetas = np.random.uniform(0, 2*np.pi, self.num_particles)
        self.directions = np.zeros((self.num_particles, 2))
        self.ratio = 10
        # ID system to assign activity value for the mixture
        #self.activityForce = np.full(self.num_particles, np.sqrt(2 * self.kB * self.temperature * self.gamma / self.tau))
        self.activityForce = np.full(self.num_particles, 5 * self.gamma)
        self.activityID = np.zeros(self.num_particles)
        if self.mixture:
            for ii in range(num_particles):
                if (ii % 2) == 0: self.activityID[ii] = 1
                else: self.activityID[ii] = 0
            self.activityForce[(self.activityID == 1)] = self.activityForce[0]*self.ratio
        # Other checks
        self.forcesContainer = []
        self.allforcesContainer = []
        self.distancesContainer = []
        self.positions_save_freq = 1000

        # Print the class instance
        print(self)

    def __str__(self):
        """Print the class variables."""
        return (
            f"{'='*50}\n"
            f"Created md object with settings:\n"
            f"Number of particles: {self.num_particles:d}\n"
            f"Temperature: {self.temperature:.1f}\n"
            f"Friction: {self.gamma:.1f}\n"
            f"Time step: {self.dt:.4f}\nBox size: Lx {self.box_size[0]:.1f} and Ly {self.box_size[1]:.1f}\n"
            f"Density: {(self.num_particles * np.pi * (self.sigma/2)**2 / (self.box_size[0]*self.box_size[1])):.1f}\n"
            f"Integrator: {self.integrator}\n"
            f"{("Potential: " + self.potentialType) if self.interaction else 'Free particles'}"
        )

    def apply_pbc(self):
        """Apply Periodic Boundary Conditions to keep particles inside the simulation box."""
        self.positions = (self.positions + self.box_size / 2) % self.box_size - self.box_size / 2

    def compute_cell_neighbours(self):
        """"Computes nearest neighbours based on cell subdivision."""

        self.neighbours = [] # Reset neighbours
        self.neighbour_counts = []
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
            self.neighbour_counts.append(len(ii_list))
        
        self.neighborCheckPositions = self.positions

    def compute_disk_neighbours(self):
        """"Computes nearest neighbours based on disk distance."""
        self.neighbour_counts[:] = 0

        for ii in range(self.num_particles):
            count = 0
            for jj in range(self.num_particles):
                if jj != ii:
                    distances = self.positions[ii] - self.positions[jj]
                    distances -= np.round(distances / self.box_size) * self.box_size
                    if np.sum(distances**2) <= (self.cutoff + self.skin)**2:
                        if count < self.max_neighbors:
                            self.neighbours[ii, count] = jj
                            count += 1
            self.neighbour_counts[ii] = count

        self.neighborCheckPositions = self.positions

    def compute_LJ_forces(self):
        """Shifted forces using LJ potential."""
        
        forceShift = (4*self.epsilon/self.cutoff) * ((12*(self.sigma/self.cutoff)**12)-(6*(self.sigma/self.cutoff)**6))
        potShift = 4 * self.epsilon * ((self.sigma/self.cutoff)**12-(self.sigma/self.cutoff)**6)
        potDerShift = (- 4 * self.epsilon * (12*((self.sigma/self.cutoff)**12)-6*((self.sigma/self.cutoff)**6))) / self.cutoff
        self.forces = np.zeros((self.num_particles, 2))  # Reset forces
        potential_energy = 0  # Reset potential
        for ii in range(self.num_particles):
            #for jj in range(ii + 1, self.num_particles):
            for k in range(self.neighbour_counts[ii]):
                jj = self.neighbours[ii, k]
                distances = self.positions[ii] - self.positions[jj]
                distances = distances - (np.round(distances/self.box_size)) * self.box_size
                distance = np.sqrt(np.sum(distances**2))
                if distance < self.cutoff:
                    potential_energy += 0.5 * (4*self.epsilon*((self.sigma/distance)**12-(self.sigma/distance)**6) - potShift - ((distance - self.cutoff)*potDerShift))
                    LJforce = ((4*self.epsilon/distance) * ((12*(self.sigma/distance)**12)-(6*(self.sigma/distance)**6))) - forceShift
                    self.forces[ii] += 0.5 * ((distances/distance) * LJforce)
                    self.forces[jj] -= 0.5 * ((distances/distance) * LJforce)
                    #self.forcesContainer.append(LJforce)
                    #self.distancesContainer.append(distance)

        self.potentialEnergy = potential_energy

    def compute_WCA_forces(self):
        """Forces using WCA potential."""
        
        self.forces = np.zeros((self.num_particles, 2))  # Reset forces
        potential_energy = np.zeros(self.num_particles)  # Reset potential
        for ii in range(self.num_particles):
            #for jj in range(ii + 1, self.num_particles):
            for k in range(self.neighbour_counts[ii]):
                jj = self.neighbours[ii, k]
                distances = self.positions[ii] - self.positions[jj]
                distances -= np.round(distances/self.box_size) * self.box_size
                distance = np.linalg.norm(distances)
                if distance < self.cutoff:
                    ratio6 = (self.sigma / distance)**6
                    ratio12 = ratio6 * ratio6
                    # 0.5 for distributing the energy in the two particles
                    potential_energy[ii] += 0.5 * 0.5 * self.epsilon * (4 * (ratio12 - ratio6) + 1)
                    potential_energy[jj] += 0.5 * 0.5 * self.epsilon * (4 * (ratio12 - ratio6) + 1)
                    WCAforce = 24 * self.epsilon * (2 * ratio12 - ratio6) / distance 
                    self.forces[ii] += 0.5 * WCAforce * distances / distance
                    self.forces[jj] -= 0.5 * WCAforce * distances / distance
                #if (ii==2) and (jj==5):
                #    self.allforcesContainer.append(WCAforce)
                #if abs(distance - self.cutoff) < 0.01:
                #    self.forcesContainer.append(WCAforce)
                #    self.distancesContainer.append(distance)
        self.potentialEnergy = np.sum(potential_energy)

    def langevin_force(self):
        """Compute stochastic White Noise and friction forces."""
        noise = np.sqrt(2 * self.kB * self.temperature * self.gamma / self.dt) * np.random.randn(self.num_particles, 2)
        return -self.gamma * self.velocities + noise

    def langevin_active_noise(self):
        """Compute Colored Noise and friction forces."""
        self.thetas += np.sqrt(2 * self.dt / self.tau) * np.random.randn(self.num_particles)
        self.directions = np.stack([np.cos(self.thetas), np.sin(self.thetas)], axis=1)
        noise = self.directions * self.self.activityForce
        return -self.gamma * self.velocities + noise

    def velocity_verlet_nve(self):
        """Velocity Verlet integration for NVE dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        self.unwrappedPositions += self.velocities * self.dt
        self.apply_pbc()
        if self.interaction:
            if self.potentialType == "WCA":
                self.compute_WCA_forces()
            elif self.potentialType == "WCAnumba":
                self.forces, self.energy = compute_WCA_forces_numba(self.positions, self.neighbours, self.neighbour_counts, 
                                                                         self.box_size, self.sigma, self.epsilon, self.cutoff)
            elif self.potentialType == "LJ":
                self.compute_LJ_forces()
        else : self.forces = np.zeros((self.num_particles, 2))
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def velocity_verlet_langevin(self):
        """Velocity Verlet integration for Langevin dynamics."""
        self.velocities += 0.5 * self.forces / self.mass * self.dt
        self.positions += self.velocities * self.dt
        self.unwrappedPositions += self.velocities * self.dt
        self.apply_pbc()
        if self.interaction:
            if self.potentialType == "WCA":
                self.compute_WCA_forces()
            elif self.potentialType == "WCAnumba":
                self.forces, self.energy = compute_WCA_forces_numba(self.positions, self.neighbours, self.neighbour_counts, 
                                                                         self.box_size, self.sigma, self.epsilon, self.cutoff)
            elif self.potentialType == "LJ":
                self.compute_LJ_forces()
        else : self.forces = np.zeros((self.num_particles, 2))
        self.forces += self.langevin_force()
        self.velocities += 0.5 * self.forces / self.mass * self.dt

    def euler_maruyama(self):
        """Euler-Maruyama integration for Active Brownian particles with Colored Noise."""
        if self.interaction:
            if self.potentialType == "WCA":
                self.compute_WCA_forces()
            elif self.potentialType == "WCAnumba":
                self.forces, self.energy = compute_WCA_forces_numba(self.positions, self.neighbours, self.neighbour_counts, 
                                                                         self.box_size, self.sigma, self.epsilon, self.cutoff)
            elif self.potentialType == "LJ":
                self.compute_LJ_forces()   
        self.thetas += np.sqrt(2 * self.dt / self.tau) * np.random.randn(self.num_particles)
        self.directions = np.stack([np.cos(self.thetas), np.sin(self.thetas)], axis=1)
        #self.positions += (self.forces / self.gamma) * self.dt + (self.activityForce / self.gamma) * self.directions * self.dt
        self.positions += (self.forces / self.gamma) * self.dt + (self.activityForce[:, None] / self.gamma) * self.directions * self.dt
        self.apply_pbc()
        #self.unwrappedPositions += (self.forces / self.gamma) * self.dt + (self.activityForce / self.gamma) * self.directions * self.dt
        self.unwrappedPositions += (self.forces / self.gamma) * self.dt + (self.activityForce[:, None] / self.gamma) * self.directions * self.dt

    def compute_potentialenergy(self):
        """Compute the potential energy of the system."""
        if self.potentialType == "WCAnumba":
            return np.sum(self.energy)
        else:
            return self.potentialEnergy

    def compute_temperature(self):
        """Compute the temperature of the system from the kinetic energy."""
        kinetic_energy = 0.5 * self.mass * np.sum(self.velocities ** 2)
        return kinetic_energy / self.num_particles
    
    def compute_kineticenergy(self):
        """Compute the kinetic energy of the system."""
        return (0.5 * self.mass * np.sum(self.velocities ** 2))