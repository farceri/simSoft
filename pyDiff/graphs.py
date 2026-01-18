import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import curve_fit
from numba import jit, njit, prange
import time

directory = '/home/auroisflying/thesis/gitVersion/simSoft/pyDiff/test/isf'

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
    plt.show()

if __name__ == '__main__':

    isf_comparison()