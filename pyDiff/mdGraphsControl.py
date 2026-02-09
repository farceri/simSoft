import os
import time
from mdFunctions import *
home = '/home/auroisflying/thesis/simSoft/pyDiff/test'

if __name__ == '__main__':
    
    start = time.time()

    particles = np.array([50], dtype=int)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([100.0], dtype=float)
    integrators = ['langevin']
    interactions = ['WCA']
    directoryList = []
    optionsDirectory = home + os.sep + f"{'em'}{'WCA'}_N{50:d}_T{1.0:.1f}_gamma{100.0:.2f}" 
    directoryList.append(optionsDirectory)

    for ii, num_particles in enumerate(particles):
        for jj, temperature in enumerate(temperatures):
            for kk, gamma in enumerate(frictions):
                for integrator in integrators:
                    for interaction in interactions:
                        
                        optionsDirectory = home + os.sep + f"{integrator}{interaction}_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}" 
                        directoryList.append(optionsDirectory)

    directory = home + os.sep + f"{'em'}{'WCA'}_N{1000:d}_T{10.0:.1f}_gamma{0.1:.2f}" + os.sep + "iteration1"
    with open(directory + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    #lastConfiguration(md, directory)

    perc=100
    densitySquares(md, directory, num_bins=6, yDivision=5, perc=perc)
    densityBands(md, directory, xDivision=100, perc=perc)

    #msdTotality(home, directoryList, title = "MSD", outputName = "msd")
    #kValues = ssfTotality(home, directoryList, title = "SSF", outputName = "ssf", graph = True)
    #taus = isfTotality(home, directoryList, kValues, title = "ISF", outputName = "isf")
    #tauPlotter(home, r"$\tau$ with increasing $T$", "tau", temperatures, taus, "seagreen")
    #cvvTotality(home, directoryList, title = "CVV", outputName = "cvv")

    print("It took %fs" %(time.time()-start))