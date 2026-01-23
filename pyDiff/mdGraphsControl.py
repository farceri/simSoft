import os
import time
from mdFunctions import *
home = '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test'

if __name__ == '__main__':
    
    start = time.time()
    temperature = 1.0
    gamma = 1.0
    dt = 0.0001
    integrator = 'nve'

    num_particles = 10
    optionsDirectory1 = home + os.sep + f"{integrator}WCA_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}" 
    num_particles = 70
    optionsDirectory2 = home + os.sep + f"{integrator}WCA_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}" 

    directoryList = {optionsDirectory1, optionsDirectory2}

    msdTotality(home, directoryList, title = "MSD", outputName = "msd")
    kValues = ssfTotality(home, directoryList, title = "SSF", outputName = "ssf", graph = True)
    isfTotality(home, directoryList, kValues, title = "ISF", outputName = "isf")

    print("It took %fs" %(time.time()-start))

# TODOS
# Implementation of taus
# Update load_data in class
# Remove the initial assignment of positions in the class for all_