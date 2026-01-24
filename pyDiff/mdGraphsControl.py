import os
import time
from mdFunctions import *
home = '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test'

if __name__ == '__main__':
    
    start = time.time()

    particles = np.array([20, 70], dtype=int)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([1.0], dtype=float)
    dt = 0.0001
    integrator = 'nve'
    interaction = 'WCA'
    directoryList = []

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):
                
                optionsDirectory = home + os.sep + f"{integrator}{interaction}_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}" 
                directoryList.append(optionsDirectory)

    msdTotality(home, directoryList, title = "MSD with different densities", outputName = "msd")
    kValues = ssfTotality(home, directoryList, title = "SSF with different densities", outputName = "ssf", graph = True)
    taus = isfTotality(home, directoryList, kValues, title = "ISF with different densities", outputName = "isf")
    #tauPlotter(home, r"$\tau$ with increasing $T$", "tau", temperatures, taus, "seagreen")
    cvvTotality(home, directoryList, title = "CVV with different densities", outputName = "cvv")

    print("It took %fs" %(time.time()-start))