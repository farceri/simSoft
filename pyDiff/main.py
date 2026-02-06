import os
import sys
import time
import pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from mdFunctions import *
from mdClasses import MolecularDynamics
#np.random.seed(0)
# python3 main.py '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test' 1e06 1e02 30 10 0.5 1e02 0/read

if __name__ == '__main__':

    start = time.time()
    directory = sys.argv[1] 
    num_steps = int(float(sys.argv[2])) 
    num_part = int(float(sys.argv[3])) # FA: ADDED NUM PARTICLES AS AN INPUT PARAMETER
    Lx = float(sys.argv[4]) # 30
    Ly = float(sys.argv[5]) # 10
    temp = float(sys.argv[6]) # 0.5
    beta = float(sys.argv[7]) # 1e02
    save_freq = int(num_steps/100)
    print_freq = int(num_steps/10)
    density = np.pi * (0.5)**2 * num_part / (Lx * Ly)
    # FA: take density and box ratio Lx/Ly as input and then compute Lx and Ly

    # Code controls
    iterations = 1
    randomizingSteps = int(1e04) # FA: easier to read than 100000
    compute_gif = False
    read_data = sys.argv[8] # FA: added as in input
    if read_data == 'read':
        load_data = True
    else:
        load_data = False
    interaction = True
    integrator = sys.argv[9] # Options are nve, langevin and em
    potentialType = sys.argv[10] # Options are LJ, WCA and WCAnumba
    print(f"Input: {integrator} integrator and {potentialType} potential - density: {density}")

    # Control the parameters 
    particles = np.array([num_part], dtype=int)
    temperatures = np.array([temp], dtype=float)
    frictions = np.array([beta], dtype=float)

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):

                # Create the path to save the data
                if interaction : optionsDirectory = f"{integrator}WCA_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}"
                else : optionsDirectory = f"{integrator}FREE_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}"

                for iteration in range(iterations):

                    # Create md object with input settings - more settings can be added
                    md = MolecularDynamics(num_particles, temperature, gamma, potentialType, interaction = interaction, 
                                           initialConf = load_data, integrator=integrator, Lx=Lx, Ly=Ly)
                    md.positions_save_freq = save_freq

                    # Create arrays for storing energy and msd
                    temp = np.empty(0)
                    msd = np.empty(0)
                    potential = np.empty(0)
                    kinetic = np.empty(0)
                    md.compute_disk_neighbours()
                    #md.compute_cell_neighbours()
                    smallOrder()
                    if load_data == True:
                        print("Reading initial configuration")
                    else:
                        print(f"Initialization: running {randomizingSteps} {integrator} steps")
                        for step in range(randomizingSteps):
                            distances = md.positions - md.neighborCheckPositions
                            distances -= np.round(distances/md.box_size) * md.box_size
                            distance = np.linalg.norm(distances, axis=1)
                            if np.any(distance >= md.skin/2): # Update the neighbour list only when necessary
                                md.compute_disk_neighbours() 
                                #md.compute_cell_neighbours()
                            if integrator == 'nve':
                                md.velocity_verlet_nve()
                            elif integrator == 'langevin':
                                md.velocity_verlet_langevin()
                            elif integrator == 'em':
                                md.euler_maruyama()

                        epot = md.compute_potentialenergy() / num_particles
                        ekin = md.compute_kineticenergy() / num_particles
                        etot = epot + ekin
                        print(f"Energy after initialization, U: {epot}, K: {ekin}, U+K: {etot}") # FA: added energy print
                        savePath = os.path.join(directory, optionsDirectory)
                        os.makedirs(savePath, exist_ok=True)
                        np.savez(os.path.join(savePath, "initialConfiguration.npz"), positions = md.positions, velocities = md.velocities)
                        # FA: CONSIDER USING THE SAME INITIAL CONFIGURATION FOR DIFFERENT VALUES OF ACTIVITY

                    # Run integration, store and print data at given frequency
                    smallOrder()
                    md.unwrappedPositions = md.positions.copy()
                    md.allPositions.append(md.positions.copy())
                    md.allUnwrappedPositions.append(md.unwrappedPositions.copy())
                    md.allVelocities.append(md.velocities.copy())
                    for step in range(num_steps + save_freq):
                        distances = md.positions - md.neighborCheckPositions
                        distances -= np.round(distances/md.box_size) * md.box_size
                        distance = np.linalg.norm(distances, axis=1)
                        if np.any(distance >= md.skin/2): # Update the neighbour list only when necessary
                            md.compute_disk_neighbours() 
                            #md.compute_cell_neighbours()
                        if integrator == 'nve':
                            md.velocity_verlet_nve()
                        elif integrator == 'langevin':
                            md.velocity_verlet_langevin()
                        elif integrator == 'em':
                            md.euler_maruyama()
                        if step % save_freq == 0:
                            temp = np.append(temp, md.compute_temperature())
                            potential = np.append(potential, md.compute_potentialenergy()/num_particles)
                            kinetic = np.append(kinetic, md.compute_kineticenergy()/num_particles)
                        if step % md.positions_save_freq == 0:
                            md.allPositions.append(md.positions.copy()) # FA: REMOVED TO OCCUPY LESS MEMORY, IT CAN BE COMPUTED IN THE ANALYSIS
                            md.allUnwrappedPositions.append(md.unwrappedPositions.copy())
                            md.allVelocities.append(md.velocities.copy())
                        if step % print_freq == 0:
                            print(f"Step {step}, T: {temp[-1]:.4f}, E: {potential[-1]+kinetic[-1]:.7f}")

                    iterationDirectory = f"iteration{iteration+1}"
                    savePath = os.path.join(directory, optionsDirectory, iterationDirectory)
                    os.makedirs(savePath, exist_ok=True)
                    # Save the class instance
                    with open(savePath + os.sep + 'classInstance.pkl', "wb") as f:
                        pickle.dump(md, f)
                    # Save time, temperature, potential and kinetic energy
                    simTime = np.arange(0, num_steps + save_freq, save_freq) * md.dt 
                    np.savetxt(savePath + os.sep + 'evolutionData.dat', np.column_stack((simTime, temp, potential, kinetic)))
                    np.savez(os.path.join(savePath, "lastConfiguration.npz"), positions = md.positions, velocities = md.velocities)

                # Plot in a gif the particles moving
                if compute_gif: 
                    md.allPositions = np.array(md.allPositions)
                    md.allPositions = np.stack(md.allPositions, axis=-1)  
                    part_evolution(md, md.allPositions)

    smallOrder()
    print("It took %fs" %(time.time()-start))
    bigOrder()