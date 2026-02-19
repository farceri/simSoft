import os
import time
from mdFunctions import *
home = '/home/auroisflying/thesis/simSoft/pyDiff/test'

if __name__ == '__main__':
    
    start = time.time()

    particles = np.array([800], dtype=int)
    densities = np.array([0.5], dtype=float)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([10.0], dtype=float)
    integrators = ['em']
    interactions = ['WCA']
    mixture = ['False']
    directoryList = []
    #optionsDirectory = home + os.sep + f"{'em'}{'WCA'}_N{800:d}_phi{0.5:.1f}_T{1.0:.1f}_g{10.0:.2f}_mixFalse" 
    #directoryList.append(optionsDirectory)

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):
                for ll, density in enumerate(densities):
                    for integrator in integrators:
                        for interaction in interactions:
                            
                            if interaction == "WCA": optionsDirectory = home + os.sep + f"{integrator}{interaction}_N{num_particles:d}_phi{density:.1f}_T{temperature:.1f}_g{gamma:.2f}" 
                            else: optionsDirectory = home + os.sep + f"{integrator}{interaction}_N{num_particles:d}_T{temperature:.1f}_g{gamma:.2f}" 
                            directoryList.append(optionsDirectory)

    directory = home + os.sep + f"tau_20.0_{'langevin'}{'WCA'}_N{1000:d}_phi{0.5:.1f}_T{1.0:.1f}_g{10:.2f}_mixFalse" + os.sep + "iteration1"
    print(directory)
    giantIdx = cluster(directory, "gcc", start=100, stop=100, howMany=1, plot=True)
    computeTemperature(directory, cluIdxs=giantIdx)
    #configuration(directory, perc=perc, outputName=f"conf{perc}", cluIdxs=0)
    #densitySquares(directory, num_bins=10, yDivision=10, perc=perc, outputName=f"hist{perc}")
    #densityBands(directory, xDivision=40, perc=perc, outputName=f"dist{perc}")

    #msdTotality(home, directoryList, outputName = "msd_Newton_freevswca")
    #kValues = ssfTotality(home, directoryList, outputName = "ssf_Lan_em", graph = True)
    #kValues = np.array((5, 5))
    #taus = isfTotality(home, directoryList, kValues, outputName = "isf_Lan_em")
    #tauPlotter(home, outputName = "tau", temperatures = temperatures, taus = taus, color = "darkslategrey")
    #cvvTotality(home, directoryList, title = "", outputName = "cvv")

    #energy_graph(home=home, directory=directoryList[0])
    #densitySquares_cluster(directory, yDivision=6, perc=perc, outputName=f"cluster{perc}", plot=True)
    #gcc_time(directory, yDivision=6, outputName=f"gcc")

    print("It took %fs" %(time.time()-start))