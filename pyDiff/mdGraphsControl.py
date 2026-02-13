import os
import time
from mdFunctions import *
home = '/home/auroisflying/thesis/simSoft/pyDiff/test'

if __name__ == '__main__':
    
    start = time.time()

    particles = np.array([1000], dtype=int)
    densities = np.array([0.5], dtype=float)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([1.0], dtype=float)
    integrators = ['em']
    interactions = ['WCA']
    directoryList = []
    #optionsDirectory = home + os.sep + f"{'langevin'}{'WCA'}_N{num_particles:d}_phi{0.5:.1f}_T{1.0:.1f}_g{100.0:.2f}" 
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

    directory = directoryList[0] + os.sep + "iteration1"
    perc = 100
    giantIdx = cluster(directory, "gcc", start=100, stop=100, howMany=1, plot=True)
    #configuration(directory, perc=perc, outputName=f"conf{perc}", cluIdxs=0)
    #densitySquares(directory, num_bins=10, yDivision=10, perc=perc, outputName=f"hist{perc}")
    #densityBands(directory, xDivision=40, perc=perc, outputName=f"dist{perc}")

    #msdTotality(home, directoryList, outputName = "msd_Newton_freevswca")
    #kValues = ssfTotality(home, directoryList, outputName = "ssf_Lan_em", graph = True)
    #kValues = np.array((5, 5))
    #taus = isfTotality(home, directoryList, kValues, outputName = "isf_Lan_em")
    #tauPlotter(home, outputName = "tau", temperatures = temperatures, taus = taus, color = "darkslategrey")
    cvvTotality(home, directoryList, title = "", outputName = "cvv")

    #densitySquares_cluster(directory, yDivision=6, perc=perc, outputName=f"cluster{perc}", plot=True)
    #gcc_time(directory, yDivision=6, outputName=f"gcc")

    print("It took %fs" %(time.time()-start))