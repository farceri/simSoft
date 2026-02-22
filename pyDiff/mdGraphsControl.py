import os
import time
from mdFunctions import *
home = '/home/auroisflying/thesis/simSoft/pyDiff/test'

if __name__ == '__main__':
    
    start = time.time()

    particles = np.array([3000], dtype=int)
    densities = np.array([0.5], dtype=float)
    temperatures = np.array([1.0], dtype=float)
    frictions = np.array([10.0], dtype=float)
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
                                    if interaction : optionsDirectory = home + os.sep + f"{integrator}WCA_N{num_particles:d}_phi{0.5:.1f}_T{temperature:.1f}_g{gamma:.2f}"
                                    else : optionsDirectory = home + os.sep + f"{integrator}FREE_N{num_particles:d}_T{temperature:.1f}_g{gamma:.2f}"

                                directoryList.append(optionsDirectory)

    directory = directoryList[0] + os.sep + "iteration1"
    giantIdx = cluster(directory, "gcc_ap_zoomed", start=0, stop=0, howMany=1, plot=True, densityStudies=False)
    #computeTemperature(directory, cluIdxs=giantIdx, perc = 100)

    #msdTotality(home, directoryList, outputName = "msd_Newton_freevswca")
    #kValues = ssfTotality(home, directoryList, outputName = "ssf_Lan_em", graph = True)
    #kValues = np.array((5, 5))
    #taus = isfTotality(home, directoryList, kValues, outputName = "isf_Lan_em")
    #tauPlotter(home, outputName = "tau", temperatures = temperatures, taus = taus, color = "darkslategrey")
    #cvvTotality(home, directoryList, title = "", outputName = "cvvComparison")

    #energy_graph(home=home, directory=directoryList[0])
    #densitySquares_cluster(directory, yDivision=6, perc=perc, outputName=f"cluster{perc}", plot=True)
    #gcc_time(directory, yDivision=6, outputName=f"gcc")

    print("It took %fs" %(time.time()-start))