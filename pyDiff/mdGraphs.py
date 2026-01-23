import os
import pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from mdFunctions import *

directory = '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test'

def isf_comparison():

    dataNVT = np.loadtxt(directory + os.sep + 'NVTdens0.39T1.0gam2.0.dat')
    dataEM = np.loadtxt(directory + os.sep + 'EMdens0.39T1.0gam2.0.dat')

    plt.title(r"comparison", fontsize=16)
    plt.plot(dataNVT[:, 0], dataNVT[:, 1], label="NVT", color="orchid", linewidth=1, linestyle='solid')
    plt.plot(dataEM[:, 0], dataEM[:, 1], label="EM", color="seagreen", linewidth=1, linestyle='solid')
    plt.xscale("log")
    plt.xlabel(r"Simulation time, $(t-t_0)$", fontsize=14)
    plt.ylabel(r"ISF", fontsize=14)
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.tight_layout()
    plt.legend()
    plt.savefig(directory + "/isfComparison.png", transparent=False, format="png")

def msdGraph(mdObject, evolutionData):

    msd = []
    msd = np.append(msd, compute_msd(mdObject.allUnwrappedPositions))

    plt.clf()
    fit = True
    if integrator == 'nve': plt.title(r"MSD with interaction:%s and Langevin:False" %(mdObject.interaction), fontsize=16)
    else : plt.title(r"MSD with interaction:%s and Langevin:True" %(mdObject.interaction), fontsize=16)

    density = (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2)
    simTime = evolutionData[:, 0]
    eq = gamma*md.dt*md.positions_save_freq

    if mdObject.integrator == 'nve': 
        if (mdObject.interaction == False) and fit:
            popt, pcov = curve_fit(fitFunc_pow, simTime, msd) 
            plt.plot(simTime, fitFunc_pow(simTime, popt[0], popt[1], popt[2]), color="seagreen", linestyle='solid', linewidth=1) 
            plt.plot(simTime, msd, color="seagreen", linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$N$=%.0f, T=%.1f, $\propto t^{%.1f}$" %(num_particles, temperature, popt[1]))
        else:
            plt.plot(simTime, msd, color="seagreen", linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$\rho$=%.2f, T=%.1f" %(density, temperature))
    elif mdObject.integrator == 'langevin' or mdObject.integrator == 'em': 
        if (mdObject.interaction == False) and fit:
            # Ballistic regime
            popt1, pcov1 = curve_fit(fitFunc_pow, simTime[:int(1/eq)], msd[:int(1/eq)]) 
            plt.plot(simTime[:int(1/eq)], fitFunc_pow(simTime[:int(1/eq)], popt1[0], popt1[1], popt1[2]), color="seagreen", linestyle='solid', linewidth=1) 
            # Diffusive regime
            popt2, pcov2 = curve_fit(fitFunc_lin, simTime[int(6/eq):], msd[int(6/eq):]) 
            plt.plot(simTime[int(6/eq):], fitFunc_lin(simTime[int(6/eq):], popt2[0], popt2[1]), color="seagreen", linestyle='solid', linewidth=1) 
            plt.plot(simTime, msd, color="seagreen", linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f, $\propto t^{%.1f}\rightarrow\propto t$" %(num_particles, temperature, gamma, popt1[1]))
        else:
            plt.plot(simTime, msd, color="seagreen", linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f" %(density, temperature, gamma))
            
    plt.ylabel(r"MSD, $\langle |r(t)-r_0|^2 \rangle$", fontsize=14)
    if integrator == 'langevin':  plt.xlabel(r"Simulation time, $t$", fontsize=14)
    else: plt.xlabel(r"Simulation time, $t$", fontsize=14)
    plt.tight_layout()
    plt.ylim(bottom=0.1, top=1.5*np.max(msd))
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig(directory + "/msdTEMP.png", transparent=False, format="png")

def ssfGraph(mdObject):

    kMods = np.linspace((2*np.pi/mdObject.box_size[0]), (4*np.pi), 30)
    ssf_total_self, ssf_total_int, chosenk = compute_ssf(mdObject.allUnwrappedPositions, mdObject.num_particles, 0, np.shape(mdObject.allUnwrappedPositions)[0]-1, kMods)

    plt.clf()
    if integrator == 'nve': plt.title(r"SSF with interaction:%s and Langevin:False" %(mdObject.interaction), fontsize=16)
    else : plt.title(r"SSF with interaction:%s and Langevin:True" %(mdObject.interaction), fontsize=16)
    density = (num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2)
    plt.plot(kMods, ssf_total_self + ssf_total_int, color="seagreen", linewidth=1, linestyle='solid', label=r"$\rho$=%.2f, T=%.1f, $\gamma=%.1f$" %(density, temperature, gamma))
    
    plt.ylabel(r"SSF", fontsize=14)
    if integrator == 'langevin':  plt.xlabel(r"|k|", fontsize=14)
    #plt.ylim(bottom=0)
    plt.axhline(y=1, color="gray", linestyle="--")
    plt.tight_layout()
    plt.legend()
    plt.savefig(directory + "/ssfTEMP.png", transparent=False, format="png")

    return chosenk

def isfGraph(mdObject, evolutionData, chosenK):

    isf_total_self, isf_total_int = compute_isf(mdObject.allUnwrappedPositions, mdObject.num_particles, 0, np.shape(mdObject.allUnwrappedPositions)[0]-1, chosenk)

    simTime = evolutionData[:, 0]

    plt.clf()
    fit = True
    if integrator == 'nve': plt.title(r"ISF with interaction:%s and Langevin:False" %(mdObject.interaction), fontsize=16)
    else : plt.title(r"ISF with interaction:%s and Langevin:True" %(mdObject.interaction), fontsize=16)
    tvalues = simTime[0:np.shape(isf_total_self)[0]]
    density = (num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2)
    eq = gamma*md.dt*mdObject.positions_save_freq

    if integrator == 'nve': 
        if (mdObject.interaction == False) and fit:
            popt, pcov = curve_fit(fitFunc_exp, tvalues, isf_total_self + isf_total_int, 
                            maxfev=100000, p0=[1, 2, 1, 0])
            plt.plot(tvalues, fitFunc_exp(tvalues, popt[0], popt[1], popt[2], popt[3]), 
                    color="seagreen", linewidth=1, linestyle='solid')
            plt.plot(tvalues, isf_total_self + isf_total_int, 
                        color="seagreen", marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$N$=%.0f, T=%.0f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(num_particles, temperature, chosenK, popt[2], popt[1]))
        else:
            plt.plot(tvalues, isf_total_self + isf_total_int, 
                        color="seagreen", marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$\rho$=%.2f, T=%.1f, $|k|=%.1f$" %(density, temperature, chosenK))
    #elif integrator == 'langevin': 
    else:
        if (mdObject.interaction == False) and fit:
            popt1, pcov1 = curve_fit(fitFunc_exp, simTime[:int(1/eq)], isf_total_self[:int(1/eq)] + isf_total_int[:int(1/eq)], 
                            maxfev=100000, p0=[1, 2, 1, 0]) 
            plt.plot(simTime[0:int(1/eq)], fitFunc_exp(simTime[:int(1/eq)], popt1[0], popt1[1], popt1[2], popt1[3]), 
                    color="seagreen", linewidth=1, linestyle='solid')
            # Diffusive regime
            #popt2, pcov2 = curve_fit(fitFunc_exp, simTime[int(8/eq):], isf_total[ii, jj, kk, int(8/eq):, mm, 0] + isf_total[ii, jj, kk, int(8/eq):, mm, 1], 
            #                maxfev=100000, p0=[1, 1, 1, 0]) 
            #plt.plot(simTime[int(8/eq):], fitFunc_exp(simTime[int(8/eq):], popt2[0], popt2[1], popt2[2], popt2[3]), 
            #        color=colors[mm], linewidth=1, linestyle='solid', alpha = transparency)
            plt.plot(tvalues, isf_total_self + isf_total_int,
                    color="seagreen", marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(num_particles, temperature, gamma, chosenK, popt1[2], popt1[1]))
        else:
            plt.plot(tvalues, isf_total_self + isf_total_int,
                        color="seagreen", marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f, $|k|=%.1f$" %(density, temperature, gamma, chosenK))
    plt.ylabel(r"ISF", fontsize=14)
    if integrator == 'langevin':  plt.xlabel(r"Simulation time, $(t-t_0)$", fontsize=14)
    else: plt.xlabel(r"Simulation time, $(t-t_0)$", fontsize=14)
    #plt.ylim(bottom=0)
    plt.xscale("log")
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.tight_layout()
    plt.legend()
    plt.savefig(directory + "/isfTEMP.png", transparent=False, format="png")

if __name__ == '__main__':
    
    num_particles = 75
    temperature = 1.0
    gamma = 1.0
    dt = 0.0001
    integrator = 'nve'

    #optionsDirectory = f"{integrator}FREE_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}"
    optionsDirectory = f"{integrator}WCA_N{num_particles:d}_T{temperature:.1f}_gamma{gamma:.2f}"
    iterationDirectory = f"iteration{1}"
    loadPath = os.path.join(directory, optionsDirectory, iterationDirectory)
    evolutionData = np.loadtxt(loadPath + os.sep + 'evolutionData.dat')

    with open(loadPath + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    msdGraph(md, evolutionData)
    chosenk = ssfGraph(md)
    isfGraph(md, evolutionData, chosenk)

# FIRST DO A CHECK FOR THE SSF AND THE ISF WITH THE OTHER CODE. IT'S PROBABLY BETTER TO DO IT WITH A FIXED SEED
# OR WITH FORCES, BECAUSE THEY WILL CHOOSE DIFFERENT K OTHERWISE.

# Implementation of iterations
# Check correctness with original code in github (save a local version)
# Implementation of taus
# Comment all functions
# Removed unused code pieces
