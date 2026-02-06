import os
import time
from mdFunctions import *
home = '/home/auroisflying/thesis/simSoft/pyDiff/test'

if __name__ == '__main__':
    
    start = time.time()

    particles = np.array([50], dtype=int)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([10.0], dtype=float)
    integrators = ['langevin']
    interactions = ['WCA']
    directoryList = []
    #optionsDirectory = home + os.sep + f"{'em'}{'FREE'}_N{25:d}_T{1.0:.1f}_gamma{10.0:.2f}" 
    #directoryList.append(optionsDirectory)

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):
                for integrator in integrators:
                    for interaction in interactions:
                        
                        optionsDirectory = home + os.sep + f"{integrator}{interaction}_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}" 
                        directoryList.append(optionsDirectory)

    #msdTotality(home, directoryList, title = "MSD", outputName = "msd")
    kValues = ssfTotality(home, directoryList, title = "SSF", outputName = "ssf", graph = True)
    print("ssf done")
    taus = isfTotality(home, directoryList, kValues, title = "ISF", outputName = "isf")
    #tauPlotter(home, r"$\tau$ with increasing $T$", "tau", temperatures, taus, "seagreen")
    #cvvTotality(home, directoryList, title = "CVV", outputName = "cvv")

    print("It took %fs" %(time.time()-start))