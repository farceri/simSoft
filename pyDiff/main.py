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
    # Read input parameters
    directory = sys.argv[1] # Directory for input and output
    num_steps = int(float(sys.argv[2])) # Number of integration steps
    save_freq = int(num_steps/100)
    print_freq = int(num_steps/10)

    # Code controls
    iterations = 1
    randomizingSteps = 5000 # Initial steps to do to randomize positions
    compute_energy = False
    compute_msdbool = True
    compute_ssfbool = True
    compute_isfbool = True
    compute_gif = False
    store_data = True 
    load_data = False # True if there is a file with initial coordinates
    if load_data: randomizingSteps = 0
    interaction = True
    integrator = 'nve' # Options are nve, langevin and eulerMaruyama
    potentialType = 'WCA' # Options are LJ and WCA

    # Control the parameters (to have a single run, just put a single value for each)
    particles = np.array([75], dtype=int)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([1.0], dtype=float)
    msd_total = np.zeros((np.shape(particles)[0], np.shape(temperatures)[0], np.shape(frictions)[0], int(num_steps/save_freq) + 1))
    isf_total = np.zeros((np.shape(particles)[0], np.shape(temperatures)[0], np.shape(frictions)[0], int(num_steps/save_freq) + 1, 1, 2))
    kvalue = np.zeros((np.shape(particles)[0], np.shape(temperatures)[0], np.shape(frictions)[0]))
    ssf_total = np.zeros((np.shape(particles)[0], np.shape(temperatures)[0], np.shape(frictions)[0], 30, 2))

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):

                for iteration in range(iterations):

                    # Create md object with input settings - more settings can be added
                    md = MolecularDynamics(num_particles, temperature, gamma, potentialType, interaction = interaction, initialConf = load_data, integrator=integrator)
                    mods = np.linspace((2*np.pi/md.box_size[0]), (4*np.pi), 30)
                    md.kmods = np.linspace((2*np.pi/md.box_size[0]), (4*np.pi), 30)
                    print("Iteration: ", iteration+1)
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
                    md.initial_positions = md.positions.copy()
                    md.unwrappedPositions = md.positions.copy()
                    md.allPositions.append(md.positions.copy())
                    md.allUnwrappedPositions.append(md.unwrappedPositions.copy())
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
                            if compute_msdbool : msd = np.append(msd, md.compute_msd())
                            potential = np.append(potential, md.potentialEnergy)
                            kinetic = np.append(kinetic, md.compute_kineticenergy())
                        if step % md.positions_save_freq == 0:
                            md.allPositions.append(md.positions.copy())
                            md.allUnwrappedPositions.append(md.unwrappedPositions.copy())
                        if step % print_freq == 0:
                            print(f"Step {step}, T: {temp[-1]:.4f}, E: {potential[-1]+kinetic[-1]:.7f}")

                    if compute_msdbool: msd_total[ii, jj, kk, :] += msd/iterations

                    if compute_ssfbool:
                        temp1, temp2 = md.compute_ssf(inf_lim = 0, sup_lim = np.shape(md.allUnwrappedPositions)[0]-1, doCycle = True)
                        ssf_total[ii, jj, kk, :, 0] += temp1/iterations
                        ssf_total[ii, jj, kk, :, 1] += temp2/iterations

                    if compute_isfbool: 
                        md.iter = iteration
                        md.chosenk = kvalue[ii, jj, kk]
                        md.compute_ssf(inf_lim = 0, sup_lim = np.shape(md.allUnwrappedPositions)[0]-1, doCycle = False)
                        temp1, temp2 = md.compute_isf(inf_lim = 0, sup_lim = np.shape(md.allUnwrappedPositions)[0]-1)
                        isf_total[ii, jj, kk, :, :, 0] += temp1/iterations
                        isf_total[ii, jj, kk, :, :, 1] += temp2/iterations
                        if iteration == 0 : kvalue[ii, jj, kk] = md.chosenk

                    if interaction : optionsDirectory = f"{integrator}WCA_N{md.num_particles:d}_T{md.temperature:.1f}_gamma{md.gamma:.2f}"
                    else : optionsDirectory = f"{integrator}FREE_N{md.num_particles:d}_T{md.temperature:.1f}_gamma{md.gamma:.2f}"
                    iterationDirectory = f"iteration{iteration+1}"
                    savePath = os.path.join(directory, optionsDirectory, iterationDirectory)
                    os.makedirs(savePath, exist_ok=True)
                    # Save the class instance
                    with open(savePath + os.sep + 'classInstance.pkl', "wb") as f:
                        pickle.dump(md, f)

                    simTime = np.arange(0, num_steps + save_freq, save_freq) * md.dt 
                    if store_data:
                        # Store time, temperature and energy in a single file
                        np.savetxt(savePath + os.sep + 'evolutionData.dat', np.column_stack((simTime, temp, potential, kinetic)))

                # Plot in a gif the particles moving
                if compute_gif: 
                    md.allPositions = np.array(md.allPositions)
                    md.allPositions = np.stack(md.allPositions, axis=-1)  
                    part_evolution(md, md.allPositions)

    smallOrder()
    print("It took %fs" %(time.time()-start))
    bigOrder()

    sim_time = np.arange(0, num_steps + save_freq, save_freq) * md.dt
    
    # Plot energy versus time (valid only for last iteration for now)
    if compute_energy:
        plt.clf()
        plt.title(r"Energies for: N=%d ($\rho$=%.1f), T=%.1f" %(md.num_particles, (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2), md.temperature), fontsize=16)
        plt.plot(simTime, kinetic, color='seagreen', linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Kinetic energy $K$")
        plt.tick_params(axis='both', labelsize=14)
        plt.plot(simTime, potential, color='steelblue', linewidth=0.9, linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Potential energy $U$")
        plt.tick_params(axis='both', labelsize=14)
        plt.plot(simTime, potential+kinetic, color='orchid', linewidth=0.9, linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Total energy $E_{tot}$")
        plt.tick_params(axis='both', labelsize=14)
        plt.ylabel("Energies", fontsize=14)
        plt.xlabel(r"Simulation time, $t$", fontsize=14)
        plt.tight_layout()
        plt.legend()
        plt.savefig(directory + "/energies.png", transparent=False, format="png")

    # Plot mean squared displacement
    if compute_msdbool:
        plt.clf()
        fit = True
        if integrator == 'nve': plt.title(r"MSD with interaction:%s and Langevin:False" %(interaction), fontsize=16)
        else : plt.title(r"MSD with interaction:%s and Langevin:True" %(interaction), fontsize=16)
        cmap = plt.get_cmap("viridis")  
        n_colors = int(np.shape(particles)[0]) * int(np.shape(temperatures)[0]) * int(np.shape(frictions)[0]) 
        colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
        counter = 0
        for ii, num_particles in enumerate(particles):
            density = (num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2)
            for jj, temperature in enumerate(temperatures):
                for kk, gamma in enumerate(frictions):
                    eq = gamma*md.dt*save_freq
                    if integrator == 'nve': 
                        if (interaction == False) and fit:
                            popt, pcov = curve_fit(fitFunc_pow, sim_time, msd_total[ii, jj, kk, :]) 
                            plt.plot(sim_time, fitFunc_pow(sim_time, popt[0], popt[1], popt[2]), color=colors[counter], linestyle='solid', linewidth=1) 
                            plt.plot(sim_time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$N$=%.0f, T=%.1f, $\propto t^{%.1f}$" %(num_particles, temperature, popt[1]))
                        else:
                            plt.plot(sim_time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$\rho$=%.2f, T=%.1f" %(density, temperature))
                    elif integrator == 'langevin' or integrator == 'eulerMaruyama': 
                        if (interaction == False) and fit:
                            # Ballistic regime
                            popt1, pcov1 = curve_fit(fitFunc_pow, sim_time[:int(1/eq)], msd_total[ii, jj, kk, :int(1/eq)]) 
                            plt.plot(sim_time[:int(1/eq)], fitFunc_pow(sim_time[:int(1/eq)], popt1[0], popt1[1], popt1[2]), color=colors[counter], linestyle='solid', linewidth=1) 
                            # Diffusive regime
                            popt2, pcov2 = curve_fit(fitFunc_lin, sim_time[int(6/eq):], msd_total[ii, jj, kk, int(6/eq):]) 
                            plt.plot(sim_time[int(6/eq):], fitFunc_lin(sim_time[int(6/eq):], popt2[0], popt2[1]), color=colors[counter], linestyle='solid', linewidth=1) 
                            plt.plot(sim_time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f, $\propto t^{%.1f}\rightarrow\propto t$" %(num_particles, temperature, gamma, popt1[1]))
                        else:
                            plt.plot(sim_time, msd_total[ii, jj, kk, :], color=colors[counter], linestyle='none', marker='o', markersize='3', fillstyle='none', 
                             label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f" %(density, temperature, gamma))
                    counter += 1
        plt.ylabel(r"MSD, $\langle |r(t)-r_0|^2 \rangle$", fontsize=14)
        if integrator == 'langevin':  plt.xlabel(r"Simulation time, $t$", fontsize=14)
        else: plt.xlabel(r"Simulation time, $t$", fontsize=14)
        plt.tight_layout()
        plt.ylim(bottom=0.1, top=1.5*np.max(msd_total))
        plt.xscale("log")
        plt.yscale("log")
        plt.legend()
        plt.savefig(directory + "/msd.png", transparent=False, format="png")

    # Plot static structure factor
    if compute_ssfbool:
        plt.clf()
        fit = True
        if integrator == 'nve': plt.title(r"SSF with interaction:%s and Langevin:False" %(interaction), fontsize=16)
        else : plt.title(r"SSF with interaction:%s and Langevin:True" %(interaction), fontsize=16)
        cmap = plt.get_cmap("viridis")  
        n_colors = int(np.shape(particles)[0]) * int(np.shape(temperatures)[0]) * int(np.shape(frictions)[0]) 
        colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
        taus = np.zeros((np.shape(temperatures)[0], np.shape(md.kmods)[0]))
        tvalues = sim_time[0:np.shape(isf_total)[3]]
        for ii, num_particles in enumerate(particles):
            density = (num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2)
            for jj, temperature in enumerate(temperatures):
                for kk, gamma in enumerate(frictions):
                    eq = gamma*md.dt*save_freq
                    plt.plot(np.linspace((2*np.pi/md.box_size[0]), (4*np.pi), 30), ssf_total[ii, jj, kk, :, 0] + ssf_total[ii, jj, kk, :, 1], 
                             color=colors[int(np.shape(frictions)[0])*int(np.shape(temperatures)[0])*ii + int(np.shape(frictions)[0])*jj + kk], 
                             linewidth=1, linestyle='solid', label=r"$\rho$=%.2f, T=%.1f, $\gamma=%.1f$" %(density, temperature, gamma))
                    
        plt.ylabel(r"SSF", fontsize=14)
        if integrator == 'langevin':  plt.xlabel(r"|k|", fontsize=14)
        #plt.ylim(bottom=0)
        plt.axhline(y=1, color="gray", linestyle="--")
        plt.tight_layout()
        plt.legend()
        plt.savefig(directory + "/ssf.png", transparent=False, format="png")

    # Plot intermediate scattering function
    if compute_isfbool:
        plt.clf()
        fit = True
        if integrator == 'nve': plt.title(r"ISF with interaction:%s and Langevin:False" %(interaction), fontsize=16)
        else : plt.title(r"ISF with interaction:%s and Langevin:True" %(interaction), fontsize=16)
        cmap = plt.get_cmap("viridis")  
        n_colors = 1
        colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
        taus = np.zeros((np.shape(temperatures)[0], np.shape(md.kmods)[0]))
        transparent = int(np.shape(particles)[0]) * int(np.shape(temperatures)[0]) * int(np.shape(frictions)[0]) 
        transparency = 1/transparent
        tvalues = simTime[0:np.shape(isf_total)[3]]
        for ii, num_particles in enumerate(particles):
            density = (num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2)
            for jj, temperature in enumerate(temperatures):
                for kk, gamma in enumerate(frictions):
                    eq = gamma*md.dt*save_freq
                    for mm, kMod in enumerate(np.atleast_1d(kvalue[ii, jj, kk])):
                        if integrator == 'nve': 
                            if (interaction == False) and fit:
                                popt, pcov = curve_fit(fitFunc_exp, tvalues, isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1], 
                                                maxfev=100000, p0=[1, 2, 1, 0])
                                plt.plot(tvalues, fitFunc_exp(tvalues, popt[0], popt[1], popt[2], popt[3]), 
                                        color=colors[mm], linewidth=1, linestyle='solid', alpha = transparency)
                                plt.plot(tvalues, isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1], 
                                         color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$N$=%.0f, T=%.0f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(num_particles, temperature, kMod, popt[2], popt[1]))
                                if np.shape(temperatures)[0] > 1: taus[jj, mm] = popt[2]
                            else:
                                plt.plot(tvalues, isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1], 
                                         color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$\rho$=%.2f, T=%.1f, $|k|=%.1f$" %(density, temperature, kMod))
                        #elif integrator == 'langevin': 
                        else:
                            if (interaction == False) and fit:
                                popt1, pcov1 = curve_fit(fitFunc_exp, simTime[:int(1/eq)], isf_total[ii, jj, kk, :int(1/eq), mm, 0] + isf_total[ii, jj, kk, :int(1/eq), mm, 1], 
                                               maxfev=100000, p0=[1, 2, 1, 0]) 
                                plt.plot(simTime[0:int(1/eq)], fitFunc_exp(simTime[:int(1/eq)], popt1[0], popt1[1], popt1[2], popt1[3]), 
                                        color=colors[mm], linewidth=1, linestyle='solid', alpha = transparency)
                                # Diffusive regime
                                #popt2, pcov2 = curve_fit(fitFunc_exp, simTime[int(8/eq):], isf_total[ii, jj, kk, int(8/eq):, mm, 0] + isf_total[ii, jj, kk, int(8/eq):, mm, 1], 
                                #                maxfev=100000, p0=[1, 1, 1, 0]) 
                                #plt.plot(simTime[int(8/eq):], fitFunc_exp(simTime[int(8/eq):], popt2[0], popt2[1], popt2[2], popt2[3]), 
                                #        color=colors[mm], linewidth=1, linestyle='solid', alpha = transparency)
                                plt.plot(tvalues, isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1],
                                        color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(num_particles, temperature, gamma, kMod, popt1[2], popt1[1]))
                                if np.shape(temperatures)[0] > 1: taus[jj, mm] = popt[2]
                            else:
                                plt.plot(tvalues, isf_total[ii, jj, kk, :, mm, 0] + isf_total[ii, jj, kk, :, mm, 1],
                                         color=colors[mm], marker='o', markersize='3', linestyle='none', fillstyle='none', alpha = transparency,
                                         label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f, $|k|=%.1f$" %(density, temperature, gamma, kMod))
                    transparency += 1/transparent
        plt.ylabel(r"ISF", fontsize=14)
        if integrator == 'langevin':  plt.xlabel(r"Simulation time, $(t-t_0)$", fontsize=14)
        else: plt.xlabel(r"Simulation time, $(t-t_0)$", fontsize=14)
        #plt.ylim(bottom=0)
        plt.xscale("log")
        plt.axhline(y=0, color="gray", linestyle="--")
        plt.tight_layout()
        plt.legend()
        plt.savefig(directory + "/isf.png", transparent=False, format="png")

        if store_data:
            for ii, num_particles in enumerate(particles):
                density = (num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2)
                for jj, temperature in enumerate(temperatures):
                    for kk, gamma in enumerate(frictions):
                        if integrator == 'langevin':
                            np.savetxt(directory + os.sep + r'isf/NVTdens%.2fT%.1fgam%.1f.dat' %(density, temperature, gamma), 
                                       np.column_stack((tvalues, isf_total[ii, jj, kk, :, 0, 0] + isf_total[ii, jj, kk, :, 0, 1])))
                        elif integrator == 'em':
                            np.savetxt(directory + os.sep + r'isf/EMdens%.2fT%.1fgam%.1f.dat' %(density, temperature, gamma), 
                                       np.column_stack((tvalues, isf_total[ii, jj, kk, :, 0, 0] + isf_total[ii, jj, kk, :, 0, 1])))

        # This needs to be changed a bit
        if np.shape(temperatures)[0] > 1:
            plt.clf()
            if integrator == 'nve': plt.title(r"$\tau$ with interaction:%s, increasing $T$" %(interaction), fontsize=16)
            else : plt.title(r"\tau$ with interaction:%s, increasing $T and Langevin" %(interaction), fontsize=16)
            for mm, kMod in enumerate(np.atleast_1d(kvalue[ii, jj, kk])):
                plt.plot(temperatures, taus[:, mm], color=colors[mm], linestyle='solid', marker='o', markersize='4', fillstyle='none', label=r"$|k|=%.1f$" %(kMod))
            plt.ylabel(r"$\tau$", fontsize=14)
            plt.xlabel(r"Temperature $T$", fontsize=14)
            plt.xscale("log")
            plt.yscale("log")
            plt.tight_layout()
            plt.legend()
            plt.savefig(directory + "/tau.png", transparent=False, format="png")
    
    # Other studies
    if False:
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