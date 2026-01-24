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
# python main.py '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test' 100000

if __name__ == '__main__':

    start = time.time()
    directory = sys.argv[1] 
    num_steps = int(float(sys.argv[2])) 
    save_freq = int(num_steps/100)
    print_freq = int(num_steps/10)

    # Code controls
    iterations = 3
    randomizingSteps = 10000 
    compute_gif = False
    load_data = False 
    if load_data: randomizingSteps = 0
    interaction = True
    integrator = 'nve' # Options are nve, langevin and em
    potentialType = 'WCA' # Options are LJ and WCA

    # Control the parameters 
    particles = np.array([20, 70], dtype=int)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([1.0], dtype=float)

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):

                # Create the path to save the data
                if interaction : optionsDirectory = f"{integrator}WCA_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}"
                else : optionsDirectory = f"{integrator}FREE_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}"

                for iteration in range(iterations):

                    # Create md object with input settings - more settings can be added
                    md = MolecularDynamics(num_particles, temperature, gamma, potentialType, interaction = interaction, initialConf = load_data, integrator=integrator)
                    md.positions_save_freq = save_freq

                    # Create arrays for storing energy and msd
                    temp = np.empty(0)
                    msd = np.empty(0)
                    potential = np.empty(0)
                    kinetic = np.empty(0)
                    md.compute_disk_neighbours()
                    #md.compute_cell_neighbours()

                    smallOrder()
                    print(f"Doing {randomizingSteps} randomizing steps...")
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
                            potential = np.append(potential, md.potentialEnergy)
                            kinetic = np.append(kinetic, md.compute_kineticenergy())
                        if step % md.positions_save_freq == 0:
                            md.allPositions.append(md.positions.copy())
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