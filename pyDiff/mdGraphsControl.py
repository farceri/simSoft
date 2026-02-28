import os
import time
import matplotlib
from mdFunctions import *
home = '/home/auroisflying/thesis/simSoft/pyDiff/test'
matplotlib.use("Agg")

if __name__ == '__main__':
    
    start = time.time()

    particles = np.array([3000], dtype=int)
    densities = np.array([0.5], dtype=float)
    temperatures = np.array([1], dtype=float)
    frictions = np.array([10], dtype=float)
    taus = np.array([20.0], dtype=float)
    active = True
    mixture = True
    activity_ratio = 0
    tau_ratio = 1
    activity = 2.5
    integrators = ['em']
    interactions = ['WCA']
    directoryList = []

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):
                for ll, density in enumerate(densities):
                    for tt, tau in enumerate(taus):
                        for integrator in integrators:
                            for interaction in interactions:

                                if active: 
                                    if mixture:
                                        if interaction : optionsDirectory = home + os.sep + f"{integrator}WCA_N{num_particles:d}_phi{0.5:.1f}_tau1{tau:.1f}_tau2{tau*tau_ratio:.1f}_v01{activity:.1f}_v02{(activity_ratio*activity):.1f}_g{gamma:.2f}"
                                        else : optionsDirectory = home + os.sep + f"{integrator}FREE_N{num_particles:d}_tau1{tau:.1f}_tau2{tau*tau_ratio:.1f}_v01{activity:.1f}_v02{(activity_ratio*activity):.1f}_g{gamma:.2f}"
                                    else:
                                        if interaction : optionsDirectory = home + os.sep + f"{integrator}WCA_N{num_particles:d}_phi{0.5:.1f}_tau{tau:.1f}_v0{activity:.1f}_g{gamma:.2f}"
                                        else : optionsDirectory = home + os.sep + f"{integrator}FREE_N{num_particles:d}_tau{tau:.1f}_v0{activity:.1f}_g{gamma:.2f}"
                                else:
                                    if interaction : optionsDirectory = home + os.sep + f"{integrator}WCA_phi{density:.1f}_T{temperature:.1f}_gamma{gamma:.2f}"
                                    else : optionsDirectory = home + os.sep + f"{integrator}FREE_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}"

                                directoryList.append(optionsDirectory)

    directory = directoryList[0] + os.sep + 'iteration1'
    #perc=10
    #giantIdx, center = cluster(directory, "gcc_0.3", start=50, stop=50, howMany=1, plot=True, LCCplot=False)
    #densityBands(directory, xDivision=30, perc=perc, outputName=f"bands{int(perc):d}", adjX=-center)
    #densitySquares(directory, num_bins=25, yDivision=10, perc=perc, outputName=f"Asquares{int(perc):d}")
    #computeTemperature(directory, cluIdxs=giantIdx, perc = 100)

    save = True
    if save:
        percs = np.array([45], dtype=float)
        bins_squares = np.array([15, 20, 25], dtype=int)
        bins_velocities = np.array([13, 20, 60], dtype=int)
        bands_divisions = np.array([35, 50, 50], dtype=int)
        adding = "ap"
        outputName_squares = "squares_" + adding
        outputName_velocities = "velocities_" + adding
        outputName_bands = "bands_" + adding
        #plt.ion()
        yLabel= False
        fig_squares, axs_squares = plt.subplots(1, 3, figsize=(12, 4))
        fig_velocities, axs_velocities = plt.subplots(1, 3, figsize=(12, 4))
        fig_bands, axs_bands = plt.subplots(1, 3, figsize=(12, 4))
        #fig, axs = plt.subplots(1, 1, figsize=(4, 4), dpi=300)
        
        for i, perc in enumerate(percs):
            if i==0: yLabel=True
            else: yLabel=False
            densitySquares(directory, num_bins=bins_squares[i], yDivision=10, perc=perc, outputName=outputName_squares, ax=axs_squares[i], yLabel=yLabel)
            giantIdx, center = cluster(directory, "gcc", start=perc, stop=perc, howMany=1, plot=False, LCCplot=False)
            velocitiesDistribution(directory, perc, cluIdxs=giantIdx, bins=bins_velocities[i], outputName=outputName_velocities, ax=axs_velocities[i], yLabel=i)
            densityBands(directory, xDivision=bands_divisions[i], perc=perc, outputName=outputName_bands, adjX=-center, ax=axs_bands[i], yLabel=yLabel)
            plt.pause(0.1)

        fig_squares.tight_layout()
        fig_velocities.tight_layout()
        fig_bands.tight_layout()
        plt.legend()
        #plt.ioff()
        fig_squares.savefig(directory + f"/{outputName_squares}.png", transparent=False, format="png", bbox_inches='tight')
        fig_velocities.savefig(directory + f"/{outputName_velocities}.png", transparent=False, format="png", bbox_inches='tight')
        fig_bands.savefig(directory + f"/{outputName_bands}.png", transparent=False, format="png", bbox_inches='tight')

    #print(directoryList)
    #msdTotality(home, directoryList, title = "", outputName = "msd_Langevin_freeT")
    #kValues = ssfTotality(home, directoryList, title="", outputName = "ssf_Langevin_freevswca", graph = True)
    #kValues = np.array((5, 5))
    #taus = isfTotality(home, directoryList, kValues, title="", outputName = "isf_comparison")
    #tauPlotter(home, outputName = "tau", temperatures = temperatures, taus = taus, color = "darkslategrey")
    #cvvTotality(home, directoryList, title = "", outputName = "cvvComparison")

    #energy_graph(home=home, directory=directoryList[0])
    #densitySquares_cluster(directory, yDivision=6, perc=perc, outputName=f"cluster{perc}", plot=True)
    #gcc_time(directory, yDivision=6, outputName=f"gcc")

    print("It took %fs" %(time.time()-start))