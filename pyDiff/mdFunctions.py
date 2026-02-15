"""
mdFunctions 
----------
These are the functions that are needed to study the MolecularDynamics class.
"""
import os
import pickle
import warnings
import numpy as np
import numba as nb
import networkx as nx
from scipy import ndimage
from scipy.spatial import ConvexHull
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.spatial.distance import pdist
import matplotlib.animation as animation
from scipy.optimize import OptimizeWarning
warnings.simplefilter("ignore", OptimizeWarning)

#-------------------------------------PRINT-VISUALS-----------------------------------------

def smallOrder() -> None:
    """Print a line of dashes "-----"."""
    print(f"{'-'*50}")

def bigOrder() -> None:
    """Print a line of double dashes "====="."""
    print(f"{'='*50}")

#-----------------------------------------FITS----------------------------------------------

def fitFunc_pow(xx: float, a: float, b: float, c: float) -> float:
    """Returns the value a*(x**b) + c."""
    return a * (xx**b) + c

def fitFunc_lin(xx: float, a: float, b: float) -> float:
    """Returns the value a*x + b."""
    return a * xx + b

def fitFunc_exp(xx: float, a: float, b: float, c: float, d: float) -> float:
    """Return the value a*e^(-(x/c)**b) + d."""
    return a * np.exp(-((xx/c)**b)) + d

#--------------------------------------OBSERVABLES------------------------------------------

def compute_msd(positions: np.ndarray) -> np.ndarray:
    """
    This function computes the Mean Squared Displacement over the timesteps.

    Parameters
    ----------
    positions : np.ndarray
        Positions of the particles. The format must be [timestep, particleID, dimension].

    Returns
    ----------
    msd : np.ndarray
        Array of the MSD meaned over all particles in formt [deltaTimestep].
    """
    initialPosition = positions[0]
    msd = np.zeros(np.shape(positions)[0]-1)
    for timestep in range(np.shape(positions)[0]-1):
        displacement = positions[timestep+1] - initialPosition
        msd[timestep] = np.mean(np.sum(displacement ** 2, axis=1))
    return msd

@nb.njit(parallel=True, fastmath=True)
def compute_ssf_jit(unwrappedPositions, kMods):
    num_timesteps, num_particles, dim = unwrappedPositions.shape
    angles = np.arange(0, 2*np.pi, np.pi/4)
    num_angles = angles.size
    
    ssf_total_self = np.ones(kMods.size, dtype=np.complex128)
    ssf_total_int = np.zeros(kMods.size, dtype=np.complex128)

    for t0 in nb.prange(num_timesteps - 1):
        for ii in range(kMods.size):
            kMod = kMods[ii]
            ssf_int = 0.0 + 0.0j
            for angle_idx in range(num_angles):
                angle = angles[angle_idx]
                kVec = np.array([kMod * np.cos(angle), kMod * np.sin(angle)])
                for i in range(num_particles):
                    for j in range(num_particles):
                        if i != j:
                            delta = unwrappedPositions[t0, i, :] - unwrappedPositions[t0, j, :]
                            ssf_int += np.exp(1j * np.dot(delta, kVec))
            ssf_total_int[ii] += ssf_int / (num_particles * num_angles)

    ssf_total_int /= (num_timesteps - 1)
    ssf_total_self = ssf_total_self.real
    ssf_total_int = ssf_total_int.real

    return ssf_total_self, ssf_total_int

def compute_ssf(unwrappedPositions: np.ndarray, kMods: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """
    This function computes the Static Structure Factor over the timesteps meaned over 8 different
    orientations of k and for different t0.

    Parameters
    ----------
    unwrappedPositions : np.ndarray
        Positions of the particles. The format must be [timestep, particleID, dimension].
    kMods : np.ndarray
        The array of modulus of k to use.

    Returns
    ----------
    ssf_total_self, ssf_total_int : tuple[np.ndarray, np.ndarray]
        The self and interacting part of the Static Structure Factor with format [kModValue].
    chosenk : float
        The value od k corresponding to the max of ssf_total_self + ssf_total_int.
    """

    num_particles = np.shape(unwrappedPositions)[1]
    angles = np.arange(0, 2*np.pi, np.pi/4)
    ssf_int = np.zeros((np.shape(kMods)[0]), dtype=complex)
    mask = ~np.eye(num_particles, dtype=bool)
    ssf_total_self = np.ones((np.shape(kMods)[0]), dtype=complex)
    ssf_total_int = np.zeros((np.shape(kMods)[0]), dtype=complex)

    inf_lim = 0
    sup_lim = np.shape(unwrappedPositions)[0]-1

    # Finding the maximum with the static structure factor
    for t0 in range(inf_lim, sup_lim): 
        ssf_int = np.zeros_like(ssf_int) 
        for ii, kMod in enumerate(kMods): 
            for angle in angles: 
                kVec = np.array((kMod * np.cos(angle), kMod * np.sin(angle))) 
                delta = unwrappedPositions[t0][:, None, :] - unwrappedPositions[t0] [None, :, :] 
                delta = delta[mask].reshape(num_particles, num_particles-1, -1) 
                ssf_int[ii] += np.sum(np.exp(1j * np.tensordot(delta, kVec, axes=([2],[0]))))/(num_particles) 
            ssf_total_int[ii] += (ssf_int[ii]) / (np.shape(angles)[0])
        
    ssf_total_int /= (sup_lim - inf_lim)
    ssf_total_self = np.real(ssf_total_self)
    ssf_total_int = np.real(ssf_total_int)

    return ssf_total_self, ssf_total_int

@nb.njit(parallel=True, fastmath=True)
def compute_isf_jit(unwrappedPositions, chosenk):

    num_timesteps, num_particles, dim = unwrappedPositions.shape
    angles = np.arange(0, 2*np.pi, np.pi/4)
    num_angles = angles.size
    isf_total_self = np.zeros(num_timesteps - 1, dtype=np.complex128)
    isf_total_int = np.zeros(num_timesteps - 1, dtype=np.complex128)

    for t0 in nb.prange(num_timesteps - 1):
        for t in range(t0, num_timesteps - 1):
            isf_self = 0.0 + 0.0j
            isf_int = 0.0 + 0.0j
            kMod = chosenk
            for angle_idx in range(num_angles):
                angle = angles[angle_idx]
                kVec = np.array([kMod * np.cos(angle), kMod * np.sin(angle)])

                for i in range(num_particles):
                    delta = unwrappedPositions[t, i, :] - unwrappedPositions[t0, i, :]
                    isf_self += np.exp(1j * np.dot(delta, kVec))/ (num_particles*num_angles)

                for i in range(num_particles):
                    for j in range(num_particles):
                        if i != j:
                            delta = unwrappedPositions[t, i, :] - unwrappedPositions[t0, j, :]
                            isf_int += np.exp(1j * np.dot(delta, kVec))/ (num_particles*(num_particles-1)*num_angles)
            isf_total_self[t-t0] += (isf_self) / ((num_timesteps - 1) - (t-t0))
            isf_total_int[t-t0] += (isf_int) / ((num_timesteps - 1) - (t-t0))

    isf_total_self = isf_total_self.real 
    isf_total_int = isf_total_int.real 

    return isf_total_self, isf_total_int

def compute_isf(unwrappedPositions: np.ndarray, chosenk: float) -> tuple[np.ndarray, np.ndarray]:
    """
    This function computes the Intermediate Scattering Function over the timesteps different t0.

    Parameters
    ----------
    unwrappedPositions : np.ndarray
        Positions of the particles. The format must be [timestep, particleID, dimension].
    chosenk : float
        The value chosen previously to compute the Intermediate Scattering Function in.

    Returns
    ----------
    isf_total_self, isf_total_int : tuple[np.ndarray, np.ndarray]
        The self and interacting part of the Intermediate Scattering Function with format [deltaTimestep].
    """

    num_particles = np.shape(unwrappedPositions)[1]
    inf_lim = 0
    sup_lim = np.shape(unwrappedPositions)[0]-1

    angles = np.arange(0, 2*np.pi, np.pi/4)
    mask = ~np.eye(num_particles, dtype=bool)
    isf_self = np.zeros((1), dtype=complex)
    isf_int = np.zeros((1), dtype=complex)
    isf_total_self = np.zeros((sup_lim), dtype=complex)
    isf_total_int = np.zeros((sup_lim), dtype=complex)

    for t0 in range(inf_lim, sup_lim):
        #print("Doing", t0)
        temporary_ip = unwrappedPositions[t0]
        for t in range(t0, sup_lim):
            isf_self = np.zeros_like(isf_self)
            isf_int = np.zeros_like(isf_int)
            for ii, kMod in enumerate(np.atleast_1d(chosenk)):
                for angle in angles:
                    kVec = np.array((kMod * np.cos(angle), kMod * np.sin(angle)))
                    isf_self[ii] += np.sum(np.exp(1j * np.matmul((unwrappedPositions[t] - temporary_ip), kVec)))/(num_particles * np.shape(angles)[0])
                    delta = unwrappedPositions[t][:, None, :] - temporary_ip[None, :, :]  
                    delta = delta[mask].reshape(num_particles, num_particles-1, -1)
                    isf_int[ii] += np.sum(np.exp(1j * np.tensordot(delta, kVec, axes=([2],[0]))))/(num_particles*(num_particles-1) * np.shape(angles)[0])
                    # isf_int[ii] += np.sum(np.exp(1j * np.matmul((self.allUnwrappedPositions[t][ :, None, :] - temporary_ip[None, :, :]), kVec)))/(self.num_particles*(self.num_particles-1) * np.shape(angles)[0])
            isf_total_self[t-t0] += (isf_self) / ((sup_lim - inf_lim) - (t-t0))
            isf_total_int[t-t0] += (isf_int) / ((sup_lim - inf_lim) - (t-t0))

    if (np.any(np.abs(np.imag(isf_total_self)) > 1e-4) or np.any(np.abs(np.imag(isf_total_int))  > 1e-4)):
        print("The imaginary part is absolutely too big.")

    isf_total_self = np.real(isf_total_self)
    isf_total_int = np.real(isf_total_int)
    
    return isf_total_self, isf_total_int

def compute_cvv(allVelocities : np.ndarray) -> np.ndarray:
    """
    This function computes the Intermediate Scattering Function over the timesteps different t0.

    Parameters
    ----------
    allVelocities : np.ndarray
        Velocities of the particles. The format must be [timestep, particleID, dimension].

    Returns
    ----------
    total_autocorrelation : np.ndarray
        The self and interacting part of the Intermediate Scattering Function with format [deltaTimestep].
    """

    num_particles = np.shape(allVelocities)[1]
    inf_lim = 0
    sup_lim = np.shape(allVelocities)[0]-1

    autocorrelation = np.zeros(1)
    total_autocorrelation = np.zeros(sup_lim)

    for t0 in range(inf_lim, sup_lim):
        temporary_iv = allVelocities[t0]
        for t in range(t0, sup_lim):
            autocorrelation = np.zeros_like(autocorrelation)
            autocorrelation += np.sum(allVelocities[t] * temporary_iv) / num_particles
            total_autocorrelation[t-t0] += (autocorrelation) / ((sup_lim - inf_lim) - (t-t0))
    
    return total_autocorrelation/total_autocorrelation[0]

#------------------------------------------GIF----------------------------------------------

def part_evolution(md, positions: np.ndarray, points: int = 1, step: int = 1) -> None:
    """
    Animation of the particles.

    Parameters
    ----------
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    positions : np.ndarray
        All the saved positions in the format [particleID, dimension, frame].
    points : int
        How many steps including the current time to show on the gif. Default is 1 (only one step).
    step : int
        If one wants to reduce the positions, this says the jump between frames. 
        Default is 1 (positions not reduced).

    Returns
    ----------
    Shows a gif of the particles. Alternitively, it can save the gif.
    """

    # The commented lines are other possible options for the gif

    fig = plt.figure()
    ax = fig.add_subplot(111)
    subdivision = 10
    plt.xticks(np.arange(-md.box_size[0]/2, md.box_size[0]/2 + md.box_size[0]/subdivision, step=md.box_size[0]/subdivision))
    plt.yticks(np.arange(-md.box_size[1]/2, md.box_size[1]/2 + md.box_size[1]/subdivision, step=md.box_size[1]/subdivision))
    #plt.grid(color = 'lightgrey', linestyle = '--', linewidth = 0.5)
    scat_dic = {}
    line_dic = {}

    # Represent the correct size of the diameter
    radius = md.sigma/2
    trans = ax.transData.transform
    inv = fig.dpi_scale_trans.inverted().transform  
    x0, y0 = trans((0,0))
    x1, y1 = trans((radius, 0))
    radius_pixels = x1 - x0
    if md.interaction : size = radius_pixels**2
    else : size = 20

    for ii in range(md.num_particles):
        # If I want one highlighted particle
        #if ii==0: scat_dic["scat{0}".format(ii)] = plt.scatter(positions[ii, 0, 0], positions[ii, 1, 0], s=size,facecolor='limegreen', edgecolor="black", linewidth=0.5)
        if md.activityID[ii] == 0:
            scat_dic["scat{0}".format(ii)] = plt.scatter(positions[ii, 0, 0], positions[ii, 1, 0], s=size, facecolor='cadetblue', edgecolor="black", linewidth=0.5)
        else: 
            scat_dic["scat{0}".format(ii)] = plt.scatter(positions[ii, 0, 0], positions[ii, 1, 0], s=size, facecolor='mediumvioletred', edgecolor="black", linewidth=0.5)
        #line_dic["line{0}".format(ii)] = plt.plot(positions[ii, 0, 0], positions[ii, 1, 0])[0]   

    # Gif visual setup    
    plt.xlim([-md.box_size[0]/2, md.box_size[0]/2])
    plt.ylim([-md.box_size[1]/2, md.box_size[1]/2])
    if md.interaction : plt.title(r"N=%d, T=%.1f, $\rho$=%.2f, $\sigma$=%.1f" %(md.num_particles, md.temperature, (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1]), md.sigma))
    else : plt.title(r"N=%d, T=%.1f" %(md.num_particles, md.temperature))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal')

    # Only take some positions, always keep the last one
    updated_positions = positions[:, :, ::step]
    temp = np.reshape(positions[:, :, -1], (positions.shape[0], positions.shape[1], 1))
    updated_positions = np.concatenate((updated_positions, temp), axis = 2)

    # Show the frames update as they go
    axtext = fig.add_axes([0.44,0.96,0,0])
    axtext.axis("off")
    time = axtext.text(0.5,0.5, str(0), ha="left", va="top")
    
    def update(frame):

        dic = {}
        data = {}
        past = 0
        past = (0 if frame < points else (frame-points))
        time.set_text("Frame: "+str(frame)+"/"+str(updated_positions.shape[2]))

        for ii in range(md.num_particles):
            dic["x{0}".format(ii)] = updated_positions[ii, 0, past:frame]
            dic["y{0}".format(ii)] = updated_positions[ii, 1, past:frame]
            data["data{0}".format(ii)] = np.stack([dic["x{0}".format(ii)], dic["y{0}".format(ii)]]).T
            scat_dic["scat{0}".format(ii)].set_offsets(data["data{0}".format(ii)])
            #line_dic["line{0}".format(ii)].set_xdata(dic["x{0}".format(ii)])
            #line_dic["line{0}".format(ii)].set_ydata(dic["y{0}".format(ii)])

        return ((scat_dic["scat{0}".format(ii)]) for ii in range(md.num_particles))
        #return ((line_dic["line{0}".format(ii)]) for ii in range(md.num_particles))

    ani = animation.FuncAnimation(fig = fig, func = update, frames = updated_positions.shape[2], blit=True)
    #ani.save('test/animation.gif', writer='imagemagick', fps=30)
    plt.show()

#-----------------------------------------OTHER---------------------------------------------

def reduce_vectors(vector1: np.ndarray, vector2: np.ndarray, eps: float) -> tuple[np.ndarray, np.ndarray]:
    """
    This function reduces the vector grouping together values that are closer than eps.
    
    Parameters
    ----------
    vector1 : np.ndarray
        Vector that gives the ordering.
    vector2 : np.ndarray
        Vector that follows the ordering of vector1.
    eps : float
        Threshold so consider the values the same.


    Returns
    ----------
    reduced_v1 : np.ndarray
        Reduced vector1.
    reduced_v2 : np.ndarray
        Reduced vector2.
    """
    # v2 is the positions
    
    vector1 = np.asarray(vector1)
    vector2 = np.asarray(vector2)

    # Sort the vector so to have the same order (of vector2)
    idx = np.argsort(vector1)
    v1_sorted = vector1[idx]
    v2_sorted = vector2[idx]
    groups = []
    current_group = [0]

    # Sort together values closer than eps
    for ii in range(1, len(v1_sorted)):
        if abs(v1_sorted[ii] - v1_sorted[current_group[-1]]) <= eps:
            current_group.append(ii)
        else:
            groups.append(current_group)
            current_group = [ii]
    groups.append(current_group)

    # Mean the groups
    reduced_v1 = np.array([v1_sorted[np.array(gg)].mean() for gg in groups])
    reduced_v2 = np.array([v2_sorted[np.array(gg)].mean() for gg in groups])

    return reduced_v1, reduced_v2

#---------------------------------------MAIN-GRAPHS-----------------------------------------

def energy_graph(home: str, directory: str) -> None:
    """
    Plot the energy evolution over time.

    Parameters
    ----------
    home : str
        Directory in which to save the plot.
    directory : str
        Directory in which to get the data.

    Returns
    ----------
    A plot comparing the potential, kinetic and total energy.
    """

    plt.clf()

    subdir = "iteration1"
    loadPath = os.path.join(directory, subdir)
    evolutionData = np.loadtxt(loadPath + os.sep + 'evolutionData.dat')
    with open(loadPath + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    #plt.title(r"Energies for: N=%d ($\rho$=%.1f), T=%.1f" %(md.num_particles, (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1]), md.temperature), fontsize=16)
    plt.plot(evolutionData[:, 0], evolutionData[:, 3], color='seagreen', linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Kinetic energy $K$")
    plt.tick_params(axis='both', labelsize=14)
    plt.plot(evolutionData[:, 0], evolutionData[:, 2], color='steelblue', linewidth=0.9, linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Potential energy $U$")
    plt.tick_params(axis='both', labelsize=14)
    plt.plot(evolutionData[:, 0], evolutionData[:, 2]+evolutionData[:, 3], color='orchid', linewidth=0.9, linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Total energy $E_{tot}$")
    plt.tick_params(axis='both', labelsize=14)
    plt.ylabel("Energies", fontsize=16)
    plt.xlabel(r"Simulation time, $t$", fontsize=16)
    plt.ylim(top=1.5)
    plt.tight_layout()
    plt.legend(fontsize=13, loc="upper center", ncols=1)
    plt.savefig(home + "/energies.png", transparent=False, format="png")

    plt.clf()
    plt.plot(evolutionData[:, 0], evolutionData[:, 2]+evolutionData[:, 3], color='mediumorchid', linewidth=0.9, linestyle='solid', marker='o', markersize='1', fillstyle='none', label="Total energy $E_{tot}$")
    plt.tick_params(axis='both', labelsize=14)
    plt.ylabel(r"$E_{tot}$", fontsize=16)
    plt.xlabel(r"Simulation time, $t$", fontsize=16)
    plt.tight_layout()
    #plt.ylim(top=0.99999)
    plt.legend(fontsize=13, loc="upper left")
    plt.savefig(home + "/energiesZoom.png", transparent=False, format="png")

def msdPlotter(simTime: np.ndarray, msd: np.ndarray, md, color: tuple) -> None:
    """
    Plot the input Mean Squared Displacement with respect to time.

    Parameters
    ----------
    simTime : np.ndarray
        The time of the simulation.
    msd : np.ndarray
        MSD array which must be the same length as simTime.
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    color : tuple
        Color of the plot.

    Returns
    ----------
    A plot showing the time evolution of Mean Squared Displacement.
    """

    density = (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1])
    eq = md.gamma*md.dt*md.positions_save_freq

    if md.integrator == 'nve': 
        if md.interaction == False:
            continuoussimTime = np.linspace(0, np.max(simTime), 1000)
            popt, pcov = curve_fit(fitFunc_pow, simTime, msd) 
            plt.plot(continuoussimTime, fitFunc_pow(continuoussimTime, popt[0], popt[1], popt[2]), color=color, linestyle='solid', linewidth=1) 
            plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"T=%.2f, $\propto t^{%.1f}$" %(md.temperature, popt[1]))
        else:
            plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"T=%.1f, $\phi$=%.2f" %(md.temperature, density))
            
    elif md.integrator == 'langevin':
        if md.interaction == False:
            # Ballistic regime
            dif = 6/eq
            if int(1/eq) > 3:
                continuoussimTime = np.linspace(0, np.max(simTime[:int(1/eq)]), 100)
                popt1, pcov1 = curve_fit(fitFunc_pow, simTime[:int(1/eq)], msd[:int(1/eq)]) 
                #plt.plot(continuoussimTime, fitFunc_pow(continuoussimTime, popt1[0], popt1[1], popt1[2]), color=color, linestyle='solid', linewidth=1) 
            # Diffusive regime
            if int(dif) < 95:
                continuoussimTime = np.linspace(np.min(simTime[int(dif):]), np.max(simTime[int(dif):]), 100)
                popt2, pcov2 = curve_fit(fitFunc_lin, simTime[int(dif):], msd[int(dif):]) 
                plt.plot(continuoussimTime, fitFunc_lin(continuoussimTime, popt2[0], popt2[1]), color=color, linestyle='solid', linewidth=1) 
            
            if int(1/eq)>3 and int(dif)<95: plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f, $\propto %.2ft$" %(md.num_particles, md.temperature, md.gamma, popt2[0]))
            elif int(dif)<95: plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f, $\propto %.2ft$" %(md.num_particles, md.temperature, md.gamma, popt2[0]))
            elif int(1/eq)>3: plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f" %(md.num_particles, md.temperature, md.gamma))
            else: plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$N$=%.0f, T=%.1f, $\gamma$=%.1f, $\rightarrow\propto t$" %(md.num_particles, md.temperature, md.gamma))
        else:
            plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
                label=r"$\phi$=%.2f, T=%.1f, $\gamma$=%.1f" %(density, md.temperature, md.gamma))
            
    elif md.integrator == 'em':
        # Diffusive regime at short time before the caging
        continuoussimTime = np.linspace(0, np.max(simTime[:int(5/eq)]), 100)
        popt, pcov = curve_fit(fitFunc_lin, simTime[:int(5/eq)], msd[:int(5/eq)]) 
        plt.plot(continuoussimTime, fitFunc_lin(continuoussimTime, popt[0], popt[1]), color=color, linestyle='solid', linewidth=1) 
        plt.plot(simTime, msd, color=color, linestyle='none', marker='o', markersize='3', fillstyle='none', 
            label=r"$\phi$=%.0f, T=%.1f, $\gamma$=%.1f, $D=%.1f$" %(density, md.temperature, md.gamma, popt[0]/4))

def msdTotality(home: str, directoryList: list[str], title: str, outputName: str) -> None:
    """
    Plot which compares different plots of the Mean Squared Displacement with means over the iterations.

    Parameters
    ----------
    home : srt
        Directory where all other sub-directories are contained.
    directoryList : list[str]
        All subdirectories to consider.
    title : srt
        The title of the plot.
    outputName : str
        The name of the plot file.

    Returns
    ----------
    A plot comparing the MSD values for different md objects.
    """

    highestValue = 0
    cmap = plt.get_cmap("viridis")  
    n_colors = int(len(directoryList)) 
    colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
    colCounter = 0

    plt.clf()
    #plt.title(title, fontsize=16)
    directoryList = sorted(directoryList)

    for directory in directoryList:
        print("\n", "Doing directory: ", directory, "\n")

        msd = None
        for root, subdirectories, files in os.walk(directory):
            subdirectories.sort()

            # Mean between the iterations
            for subdir in subdirectories:
                print("Doing subdirectory: ", subdir)

                loadPath = os.path.join(directory, subdir)
                evolutionData = np.loadtxt(loadPath + os.sep + 'evolutionData.dat')
                with open(loadPath + os.sep +"classInstance.pkl", "rb") as f:
                    md = pickle.load(f)

                if msd is None: msd = compute_msd(md.allUnwrappedPositions)/np.shape(subdirectories)[0]
                else: msd += compute_msd(md.allUnwrappedPositions)/np.shape(subdirectories)[0]

        highestValue = max(highestValue, np.max(msd))
        msdPlotter(evolutionData[:, 0], msd, md, colors[colCounter])
        colCounter += 1

    plt.ylabel(r"MSD, $\langle |r(t)-r_0|^2 \rangle$", fontsize=14)
    plt.xlabel(r"Simulation time, $t$", fontsize=14)
    plt.tight_layout()
    plt.ylim(bottom=0.01, top=1.5*highestValue)
    plt.xlim(left=0.05, right=12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig(home + f"/{outputName}.png", transparent=False, format="png")

def ssfPlotter(kMods: np.ndarray, ssf_self: np.ndarray, ssf_int: np.ndarray, md, color: tuple) -> None:
    """
    Plot the input Static Structure Factor with respect to magnitudes of k.

    Parameters
    ----------
    kMods : np.ndarray
        The magnitudes of k values to use.
    ssf_self : np.ndarray
        SSF array of the self part.
    ssf_int : np.ndarray
        SSF array of the interacting part.
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    color : tuple
        Color of the plot.

    Returns
    ----------
    A plot showing the Static Structure Factor with increasing k magnitudes.
    """

    density = (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1])
    if md.integrator == 'nve': 
        if md.interaction == False:
            plt.plot(kMods, ssf_self + ssf_int, color=color, linewidth=1, linestyle='solid', label=r"$N$=%.0f, T=%.1f" %(md.num_particles, md.temperature))
        else:
            plt.plot(kMods, ssf_self + ssf_int, color=color, linewidth=1, linestyle='solid', label=r"$\rho$=%.2f, T=%.1f" %(density, md.temperature))
    else:
        if md.interaction == False:
            plt.plot(kMods, ssf_self + ssf_int, color=color, linewidth=1, linestyle='solid', label=r"$N$=%.0f, T=%.1f, $\gamma=%.1f$" %(md.num_particles, md.temperature, md.gamma))
        else:
            plt.plot(kMods, ssf_self + ssf_int, color=color, linewidth=1, linestyle='solid', label=r"$\rho$=%.2f, T=%.1f, $\gamma=%.1f$" %(density, md.temperature, md.gamma))
    
def ssfTotality(home: str, directoryList: list[str], title: str, outputName: str, graph: bool) -> np.ndarray:
    """
    Plot which compares different plots of the Static Structure Factor with means over the iterations.

    Parameters
    ----------
    home : srt
        Directory where all other sub-directories are contained.
    directoryList : list[str]
        All subdirectories to consider.
    title : srt
        The title of the plot.
    outputName : str
        The name of the plot file.
    graph : bool
        If true, saves a plot.

    Returns
    ----------
    kValues : np.ndarray
        An array containing the values of k for which SSF is max.
    """

    highestValue = 0
    dirCounter = 0
    kValues = np.zeros(len(directoryList))
    cmap = plt.get_cmap("viridis")  
    n_colors = int(len(directoryList)) 
    colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
    colCounter = 0

    plt.clf()
    #if graph : plt.title(title, fontsize=16)
    directoryList = sorted(directoryList)

    for directory in directoryList:
        print("Doing directory: ", directory)
        
        ssf_self = None
        ssf_int = None
        for root, subdirectories, files in os.walk(directory):
            subdirectories.sort()

            # Mean between the iterations
            for subdir in subdirectories:
                print("Doing subdirectory: ", subdir)

                loadPath = os.path.join(directory, subdir)
                evolutionData = np.loadtxt(loadPath + os.sep + 'evolutionData.dat')
                with open(loadPath + os.sep +"classInstance.pkl", "rb") as f:
                    md = pickle.load(f)

                maxSize = max(md.box_size[0], md.box_size[1])
                kMods = np.linspace((2*np.pi/maxSize), (4*np.pi), 100)
                unwrappedPositions = np.array(md.allUnwrappedPositions, dtype=np.float64)
                ssf_self_temp, ssf_int_temp = compute_ssf_jit(unwrappedPositions, kMods)
                if ssf_self is None: 
                    ssf_self = ssf_self_temp/np.shape(subdirectories)[0]
                    ssf_int = ssf_int_temp/np.shape(subdirectories)[0]
                else: 
                    ssf_self += ssf_self_temp/np.shape(subdirectories)[0]
                    ssf_int += ssf_int_temp/np.shape(subdirectories)[0]

        highestValue = max(highestValue, np.max(ssf_self + ssf_int))
        chosenk = np.array([kMods[np.argmax(ssf_self + ssf_int)]])
        kValues[dirCounter] = chosenk[0]
        dirCounter += 1
        if graph : ssfPlotter(kMods, ssf_self, ssf_int, md, colors[colCounter])
        colCounter += 1

    if graph: 
        plt.ylabel(r"SSF", fontsize=14)
        plt.xlabel(r"|k|", fontsize=14)
        #plt.ylim(bottom=0)
        plt.axhline(y=1, color="gray", linestyle="--")
        plt.tight_layout()
        plt.legend()
        plt.savefig(home + f"/{outputName}.png", transparent=False, format="png")

    return kValues

def isfPlotter(simTime: np.ndarray, isf_self: np.ndarray, isf_int: np.ndarray, md, chosenk: float, color: tuple) -> float:
    """
    Plot the input Intermediate Scattering Function with respect to time.

    Parameters
    ----------
    kMods : np.ndarray
        The magnitudes of k values to use.
    isf_self : np.ndarray
        ISF array of the self part.
    isf_int : np.ndarray
        ISF array of the interacting part.
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    chosenk : float
        The value the ISF is computed in.
    color : tuple
        Color of the plot.

    Returns
    ----------
    tau : float
        Value corresponding to the ISF fit.
    """

    density = (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1])
    tau = 0

    mask = (isf_self + isf_int) > 0.2
    continuoussimTime = np.linspace(0, np.max(simTime[mask]), 10000)
    popt, pcov = curve_fit(fitFunc_exp, simTime[mask], isf_self[mask] + isf_int[mask], maxfev=10000000, p0=[1, 2, 1, 0], bounds=([0, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf]))
    plt.plot(continuoussimTime, fitFunc_exp(continuoussimTime, popt[0], popt[1], popt[2], popt[3]), color=color, linewidth=1, linestyle='solid')
    tau = popt[2]

    if md.integrator == 'nve': 
        if md.interaction == False:
            plt.plot(simTime, isf_self + isf_int, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$N$=%.0f, T=%.0f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(md.num_particles, md.temperature, chosenk, popt[2], popt[1]))
        else:
            plt.plot(simTime, isf_self + isf_int, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$\rho$=%.2f, T=%.1f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(density, md.temperature, chosenk, popt[2], popt[1]))
    else:
        if md.interaction == False:
            plt.plot(simTime, isf_self + isf_int, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$N$=%.0f, T=%.2f, $\gamma$=%.1f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(md.num_particles, md.temperature, md.gamma, chosenk, popt[2], popt[1]))
        else:
            plt.plot(simTime, isf_self + isf_int, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f, $|k|=%.1f$, $\propto e^{-((t-t_0)/%.1f)^{%.1f}}$" %(density, md.temperature, md.gamma, chosenk, popt[2], popt[1]))
    
    return tau

def isfPlotter_sep(simTime: np.ndarray, isf_self: np.ndarray, isf_int: np.ndarray, md, chosenk: float, color: tuple) -> float:
    """
    Plot the input Intermediate Scattering Function with respect to time.

    Parameters
    ----------
    kMods : np.ndarray
        The magnitudes of k values to use.
    isf_self : np.ndarray
        ISF array of the self part.
    isf_int : np.ndarray
        ISF array of the interacting part.
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    chosenk : float
        The value the ISF is computed in.
    color : tuple
        Color of the plot.

    Returns
    ----------
    tau : float
        Value corresponding to the ISF fit.
    """

    density = (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1])
    tau = 0

    mask = (isf_self + isf_int) > 0.2
    continuoussimTime = np.linspace(0, np.max(simTime[mask]), 10000)
    popt, pcov = curve_fit(fitFunc_exp, simTime[mask], isf_self[mask] + isf_int[mask], maxfev=10000000, p0=[1, 2, 1, 0], bounds=([0, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf]))
    #plt.plot(continuoussimTime, fitFunc_exp(continuoussimTime, popt[0], popt[1], popt[2], popt[3]), color=color, linewidth=1, linestyle='solid')
    tau = popt[2]

    if md.integrator == 'nve': 
        if md.interaction == False:
            plt.plot(simTime, isf_self, color=color, linestyle='solid', fillstyle='none', label=r"$I_{s}$")
            plt.plot(simTime, isf_int, color=color, linestyle='dashed', fillstyle='none',label=r"$I_{int}$")
        else:
            plt.plot(simTime, isf_self, color=color, linestyle='solid', fillstyle='none',label=r"$I_{s}$")
            plt.plot(simTime, isf_int, color=color, linestyle='dashed', fillstyle='none', label=r"$I_{int}$")
    else:
        if md.interaction == False:
            plt.plot(simTime, isf_self, color=color, linestyle='solid', fillstyle='none', label=r"$I_{s}$")
            plt.plot(simTime, isf_int, color=color, linestyle='dashed', fillstyle='none', label=r"$I_{int}$")
        else:
            plt.plot(simTime, isf_self, color=color, linestyle='solid', fillstyle='none', label=r"$I_{s}$")
            plt.plot(simTime, isf_int, color=color,  linestyle='dashed', fillstyle='none', label=r"$I_{int}$")

    return tau
   
def isfTotality(home: str, directoryList: list[str], kValues: np.ndarray, title: str, outputName: str) -> np.ndarray:
    """
    Plot which compares different plots of the Intermediate Scattering Function with means over the iterations.

    Parameters
    ----------
    home : srt
        Directory where all other sub-directories are contained.
    directoryList : list[str]
        All subdirectories to consider.
    kValues : np.ndarray
        The k magnitudes at which to compute the ISF.
    title : srt
        The title of the plot.
    outputName : str
        The name of the plot file.

    Returns
    ----------
    taus : np.ndarray
        Array containing the fit values for tau for each directory.
    """

    highestValue = 0
    dirCounter = 0
    cmap = plt.get_cmap("viridis")  
    n_colors = int(len(directoryList)) 
    colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
    colCounter = 0
    taus = np.zeros(len(directoryList))

    plt.clf()
    #plt.title(title, fontsize=16)
    directoryList = sorted(directoryList)

    for directory in directoryList:
        print("Doing directory: ", directory)
        
        isf_self = None
        isf_int = None
        for root, subdirectories, files in os.walk(directory):
            subdirectories.sort()

            # Mean between the iterations
            for subdir in subdirectories:
                print("Doing subdirectory: ", subdir)

                loadPath = os.path.join(directory, subdir)
                evolutionData = np.loadtxt(loadPath + os.sep + 'evolutionData.dat')
                with open(loadPath + os.sep +"classInstance.pkl", "rb") as f:
                    md = pickle.load(f)

                unwrappedPositions = np.array(md.allUnwrappedPositions, dtype=np.float64)
                isf_self_temp, isf_int_temp = compute_isf_jit(unwrappedPositions, kValues[dirCounter])
                if isf_self is None: 
                    isf_self = isf_self_temp/np.shape(subdirectories)[0]
                    isf_int = isf_int_temp/np.shape(subdirectories)[0]
                else: 
                    isf_self += isf_self_temp/np.shape(subdirectories)[0]
                    isf_int += isf_int_temp/np.shape(subdirectories)[0]

        highestValue = max(highestValue, np.max(isf_self + isf_int))
        #taus[dirCounter] = isfPlotter(evolutionData[:, 0], isf_self, isf_int, md, kValues[dirCounter], colors[colCounter])
        taus[dirCounter] = isfPlotter(evolutionData[:, 0], isf_self, isf_int, md, np.mean(kValues), colors[colCounter])
        colCounter += 1
        dirCounter += 1

    plt.ylabel(r"ISF", fontsize=14)
    plt.xlabel(r"Simulation time, $(t-t_0)$", fontsize=14)
    plt.xlim(left=0.02)
    plt.ylim(top=1.4)
    plt.xscale("log")
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.tight_layout()
    plt.legend()
    plt.savefig(home + f"/{outputName}.png", transparent=False, format="png")

    return taus

def tauPlotter(home: str, title: str, outputName: str, temperatures: np.ndarray, taus: np.ndarray, color: tuple) -> None:
    """
    Plot which compares different plots of the Intermediate Scattering Function with means over the iterations.

    Parameters
    ----------
    home : srt
        Directory where all other sub-directories are contained.
    title : srt
        The title of the plot.
    outputName : str
        The name of the plot file.
    temperatures : np.ndarray
        Temperatures to put on the x axis.
    taus : np.ndarray
        Taus to put on the y axis.
    color: tuple
        Color of the graph.

    Returns
    ----------
    A plot with tau changing over the temperatures.
    """

    plt.clf()
    #plt.title(title, fontsize=16)
    plt.plot(temperatures, taus, color=color, linestyle='solid', marker='o', markersize='4', fillstyle='none')
    plt.ylabel(r"$\tau$", fontsize=14)
    plt.xlabel(r"Temperature $T$", fontsize=14)
    plt.xscale("log")
    plt.yscale("log")
    plt.tight_layout()
    #plt.legend()
    plt.savefig(home + f"/{outputName}.png", transparent=False, format="png")

def cvvPlotter(simTime: np.ndarray, cvv: np.ndarray, md, color: tuple) -> float:
    """
    Plot the input velocity autocorrelation function with respect to time.

    Parameters
    ----------
    kMods : np.ndarray
        The magnitudes of k values to use.
    cvv : np.ndarray
        Velocity auto-correlation function.
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    color : tuple
        Color of the plot.

    Returns
    ----------
    tau : float
        Value corresponding to the cvv fit.
    """

    density = (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1])
    mask = (cvv) > 0.2
    continuoussimTime = np.linspace(0, np.max(simTime[mask]), 10000)

    popt, pcov = curve_fit(fitFunc_exp, simTime[mask], cvv[mask], maxfev=100000, p0=[1, 2, 1, 0])
    plt.plot(continuoussimTime, fitFunc_exp(continuoussimTime, popt[0], popt[1], popt[2], popt[3]), color=color, linewidth=1, linestyle='solid')
    tau = popt[2]

    if md.integrator == 'nve': 
        if md.interaction == False:
            plt.plot(simTime, cvv, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$N$=%.0f, T=%.0f, $\propto e^{-((t-t_0)/%.2f)^{%.3f}}$" %(md.num_particles, md.temperature, popt[2], popt[1]))
        else:
            plt.plot(simTime, cvv, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$\rho$=%.2f, T=%.1f, $\propto e^{-((t-t_0)/%.2f)^{%.3f}}$" %(density, md.temperature, popt[2], popt[1]))
    else:
        if md.interaction == False:
            plt.plot(simTime, cvv, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$N$=%.0f, T=%.0f, $\gamma$=%.1f, $\propto e^{-((t-t_0)/%.3f)^{%.1f}}$" %(md.num_particles, md.temperature, md.gamma, popt[2], popt[1]))
        else:
            plt.plot(simTime, cvv, color=color, marker='o', markersize='3', linestyle='none', fillstyle='none', 
                        label=r"$\rho$=%.2f, T=%.1f, $\gamma$=%.1f, $\propto e^{-((t-t_0)/%.3f)^{%.1f}}$" %(density, md.temperature, md.gamma, popt[2], popt[1]))

    return tau

def cvvTotality(home: str, directoryList: list[str], title: str, outputName: str) -> np.ndarray:
    """
    Plot which compares different plots of the Intermediate Scattering Function with means over the iterations.

    Parameters
    ----------
    home : srt
        Directory where all other sub-directories are contained.
    directoryList : list[str]
        All subdirectories to consider.
    title : srt
        The title of the plot.
    outputName : str
        The name of the plot file.

    Returns
    ----------
    taus : np.ndarray
        Array containing the fit values for tau for each directory.
    """

    highestValue = 0
    dirCounter = 0
    cmap = plt.get_cmap("viridis")  
    n_colors = int(len(directoryList)) 
    colors = [cmap(ii / (n_colors)) for ii in range(n_colors)]
    colCounter = 0
    taus = np.zeros(len(directoryList))

    plt.clf()
    #plt.title(title, fontsize=16)
    directoryList = sorted(directoryList)

    for directory in directoryList:
        #print("\n", "Doing directory: ", directory, "\n")
        
        cvv = None
        for root, subdirectories, files in os.walk(directory):
            subdirectories.sort()

            # Mean between the iterations
            for subdir in subdirectories:
                #print("Doing subdirectory: ", subdir)

                loadPath = os.path.join(directory, subdir)
                evolutionData = np.loadtxt(loadPath + os.sep + 'evolutionData.dat')
                with open(loadPath + os.sep +"classInstance.pkl", "rb") as f:
                    md = pickle.load(f)

                if cvv is None: cvv = compute_cvv(md.allVelocities)/np.shape(subdirectories)[0]
                else: cvv += compute_cvv(md.allVelocities)/np.shape(subdirectories)[0]

        highestValue = max(highestValue, np.max(cvv))
        taus[dirCounter] = cvvPlotter(evolutionData[:, 0], cvv, md, colors[colCounter])
        colCounter += 1
        dirCounter += 1

    plt.ylabel(r"$C_{vv}$", fontsize=14)
    plt.xlabel(r"Simulation time, $(t-t_0)$", fontsize=14)
    plt.xlim(left=0.01)
    #plt.ylim(bottom=0)
    plt.xscale("log")
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.tight_layout()
    plt.legend()
    plt.savefig(home + f"/{outputName}.png", transparent=False, format="png")

def configuration(directory, perc, outputName, cluIdxs = 0, exIdxs=0, adjX=0):

    with open(directory + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    step = int(np.shape(md.allPositions)[0]*perc/100)
    if step >= np.shape(md.allPositions)[0]: step = np.shape(md.allPositions)[0]-1
    positions = md.allPositions[step]

    if cluIdxs == 0: cluIdxs = np.arange(0, md.num_particles)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    subdivision = 5
    plt.xticks(np.arange(-md.box_size[0]/2, md.box_size[0]/2 + md.box_size[0]/subdivision, step=md.box_size[0]/subdivision))
    plt.yticks(np.arange(-md.box_size[1]/2, md.box_size[1]/2 + md.box_size[1]/subdivision, step=md.box_size[1]/subdivision))

    #data = np.load(directory + os.sep + 'lastConfiguration.npz')
    #positions = data["positions"]
    adjY = 0

    plt.xlim([-md.box_size[0]/2, md.box_size[0]/2])
    plt.ylim([-md.box_size[1]/2, md.box_size[1]/2])

    if md.interaction : plt.title(r"T=%.1f, $\rho$=%.2f, $\gamma$=%.1f, $t$=%.1f" %(md.temperature, (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]*md.box_size[1]), md.gamma, int(md.steps*md.dt*perc/100)))
    else : plt.title(r"N=%d, T=%.1f" %(md.num_particles, md.temperature))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal')

    radius = md.sigma/2
    trans = ax.transData.transform
    inv = fig.dpi_scale_trans.inverted().transform  
    x0, y0 = trans((0,0))
    x1, y1 = trans((radius, 0))
    radius_pixels = x1 - x0
    if md.interaction : size = radius_pixels**2
    else : size = 20

    positions[:, 0] += adjX
    positions[:, 0] = (positions[:, 0] + md.box_size[0] / 2) % md.box_size[0] - md.box_size[0] / 2
    positions[:, 1] += adjY
    positions[:, 1] = (positions[:, 1] + md.box_size[1] / 2) % md.box_size[1] - md.box_size[1] / 2

    for ii in range(md.num_particles):
        if md.activityID[ii] == 0:
            if ii not in cluIdxs:
                plt.scatter(positions[ii, 0], positions[ii, 1], s=size, facecolor='lightblue', edgecolor="skyblue", linewidth=0.5)
            else:
                plt.scatter(positions[ii, 0], positions[ii, 1], s=size, facecolor='cadetblue', edgecolor="black", linewidth=0.5)
        if md.activityID[ii] == 1:
            if ii not in cluIdxs:
                plt.scatter(positions[ii, 0], positions[ii, 1], s=size, facecolor='thistle', edgecolor="orchid", linewidth=0.5)
            else:
                plt.scatter(positions[ii, 0], positions[ii, 1], s=size, facecolor='mediumvioletred', edgecolor="black", linewidth=0.5)

    plt.savefig(directory + f"/{outputName}.png", transparent=False, format="png")

def densitySquares(directory, num_bins, yDivision, perc, outputName):

    with open(directory + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    step = int(np.shape(md.allPositions)[0]*perc/100)
    if step >= np.shape(md.allPositions)[0]: step = np.shape(md.allPositions)[0]-1
    positions = md.allPositions[step]

    sideSquares = md.box_size[1]/yDivision

    x = -md.box_size[0]/2
    y = -md.box_size[1]/2

    density = np.zeros((yDivision * int(md.box_size[0] / md.box_size[1]), yDivision), dtype=int)

    for ii in range(yDivision):
        x = -md.box_size[0]/2
        for jj in range(yDivision * int(md.box_size[0] / md.box_size[1])):
            for particle in range(md.num_particles):
                if (x <= positions[particle, 0] < x + sideSquares) and (y <= positions[particle, 1] < y + sideSquares):
                    density[jj, ii] += 1
            x += sideSquares
        y += sideSquares

    density = density*(np.pi*(0.5)**2)/sideSquares**2
    dens = density.ravel() 

    hist, bin_edges = np.histogram(dens, bins=num_bins, density=True)  
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    plt.figure()
    plt.bar(bin_centers, hist, width=bin_edges[1]-bin_edges[0], align='center', color="darkslategrey")
    plt.xlabel(r"$\phi$")
    plt.ylabel(r"$P(\phi)$")
    plt.tight_layout()
    #plt.legend()
    plt.savefig(directory + f"/{outputName}.png", transparent=False, format="png")

def densityBands(directory, xDivision, perc, outputName):

    with open(directory + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)
    adjX = 26

    step = int(np.shape(md.allPositions)[0]*perc/100)
    if step >= np.shape(md.allPositions)[0]: step = np.shape(md.allPositions)[0]-1
    positions = md.allPositions[step]

    sideBands = md.box_size[0]/xDivision

    x = -md.box_size[0]/2

    density = np.zeros((xDivision), dtype=int)
    positions[:, 0] += adjX
    positions[:, 0] = (positions[:, 0] + md.box_size[0] / 2) % md.box_size[0] - md.box_size[0] / 2

    for ii in range(xDivision):
        for particle in range(md.num_particles):
            if (x <= positions[particle, 0] < x + sideBands):
                density[ii] += 1
        x += sideBands

    density = density*(np.pi*(0.5)**2)/(sideBands*md.box_size[1])

    plt.figure()
    plt.plot(np.linspace(-md.box_size[0]/2, md.box_size[0]/2, xDivision), density, color="darkslategrey")
    plt.xlabel(r"$L_x$")
    plt.ylabel(r"$\phi$")
    plt.tight_layout()
    #plt.legend()
    plt.savefig(directory + f"/{outputName}.png", transparent=False, format="png")

def cluster(directory, outputName, start, stop, howMany, plot):

    percs = np.linspace(start, stop, howMany)
    gccsize = np.zeros(np.shape(percs)[0])

    with open(directory + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    threshold = 1.1
    graph = nx.Graph()
    graph.add_nodes_from(range(md.num_particles))

    for pp, perc in enumerate(percs):

        print("Perc: ", perc)
        graph.clear_edges()
        step = int(np.shape(md.allPositions)[0]*perc/100)
        if step >= np.shape(md.allPositions)[0]: step = np.shape(md.allPositions)[0]-1
        positions = md.allPositions[step]

        for ii in range(md.num_particles):
            for jj in range(ii+1, md.num_particles):
                distances = positions[ii] - positions[jj]
                distances -= np.round(distances/md.box_size) * md.box_size
                distance = np.linalg.norm(distances)
                if distance < threshold:
                    graph.add_edge(ii, jj)

        components = list(nx.connected_components(graph))
        giant_size = max(len(cc) for cc in components) if components else 0
        gccsize[pp] = giant_size / md.num_particles
        giant_component = max(components, key=len)
        giant_indices = sorted(giant_component)

        """external = []
        threshold2 = 1.5
        for ii in giant_indices:
            for jj in np.setdiff1d(range(md.num_particles), giant_indices):
                distances = positions[ii] - positions[jj]
                distances -= np.round(distances/md.box_size) * md.box_size
                distance = np.linalg.norm(distances)
                if distance < threshold2:
                    external.append(ii)"""

        # Finding the center of the cluster:
        sum1 = 0
        sum2 = 0
        for ii in giant_indices:
            theta = 2 * np.pi * ((positions[ii, 0] + md.box_size[0]/2) % md.box_size[0]) / md.box_size[0]
            sum1 += np.sin(theta)
            sum2 += np.cos(theta)

        center = md.box_size[0] * np.arctan2(sum1, sum2) / (2 * np.pi)
        center = (center + md.box_size[0]/2) % md.box_size[0] - md.box_size[0]/2

        if plot: configuration(directory, perc=perc, outputName=f"conf{int(perc):d}", cluIdxs=giant_indices, adjX=-center)

    plt.figure()
    plt.plot(percs, gccsize, color="darkslategrey")
    plt.xlabel(r"perc")
    plt.ylabel(r"size GCC")
    plt.tight_layout()
    #plt.legend()
    plt.savefig(directory + f"/{outputName}.png", transparent=False, format="png")

    return giant_indices

#--------------------------------------OTHER-GRAPHS-----------------------------------------

def plotPotForce(directory):

    epsilon = 1
    sigma = 1
    distw = np.linspace(0.9, 2**(1/6), 1000)
    distl = np.linspace(0.9, 1.8, 1000)
    dist3 = np.linspace(2**(1/6), 1.8, 1000)

    plt.plot(distw, 4*epsilon*((sigma/distw)**12-(sigma/distw)**6) + epsilon, label="WCA potential", color="mediumorchid")
    plt.plot(dist3, 0*dist3, color="mediumorchid")
    plt.plot(distl, 4*epsilon*((sigma/distl)**12-(sigma/distl)**6), label="LJ potential", color="seagreen")
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.xlabel(r"$r_{ij}$", fontsize=16)
    plt.ylabel(r"$U_{WCA}$", fontsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.tick_params(axis='both', labelsize=14)
    plt.savefig(directory + "/WCApotential.png", transparent=False, format="png")

    plt.clf()
    plt.plot(distw, ((4*epsilon/distw) * ((12*(sigma/distw)**12)-(6*(sigma/distw)**6))), label="WCA force", color="orchid")
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.xlabel(r"$r_{ij}$")
    plt.ylabel(r"$F_{WCA}$")
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.tick_params(axis='both', labelsize=14)
    plt.savefig(directory + "/WCAforce.png", transparent=False, format="png")


def forcesPotential(directory: str) -> None:
    """
    Plot two graphs: one above of the used potential and one below of the theoretical force
    compared to the actual forces during the simulation.

    Parameters
    ----------
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    directory : str
        Directory in which to save the plot.

    Returns
    ----------
    A plot of the potential and the forces.
    """

    with open(directory + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    fig, ax = plt.subplots(2, 1, figsize = (7, 7), sharex = True, dpi = 120)
    dist = np.linspace(md.sigma*0.9, md.cutoff, 1000)
    added = np.linspace(md.cutoff, md.cutoff*1.1, 1000)
    #forces, distances = reduce_vectors(md.forcesContainer, md.distancesContainer, 1e-06)
    ax[0].axhline(y=0, color="gray", linestyle="--")
    ax[0].axvline(x=md.cutoff, color="gray", linestyle="--")
    ax[1].axhline(y=0, color="gray", linestyle="--")
    ax[1].axvline(x=md.cutoff, color="gray", linestyle="--")
    if md.potentialType == "WCA" or md.potentialType == "WCAnumba":
        ax[0].plot(dist, 4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6) + md.epsilon, label="WCA potential", color="orchid")
        ax[0].plot(added, 0*added, color="orchid")
        ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6))), label="WCA force", color="orchid")
        ax[1].plot(added, 0*added, color="orchid")
    elif md.potentialType == "LJ":
        ax[0].plot(dist, 4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6), label="LJ potential")
        ax[0].plot(dist, (4*md.epsilon*((md.sigma/dist)**12-(md.sigma/dist)**6)) - (4*md.epsilon*((md.sigma/md.cutoff)**12-(md.sigma/md.cutoff)**6)) + (dist-md.cutoff)*((4*md.epsilon/md.cutoff) * ((12*(md.sigma/md.cutoff)**12)-(6*(md.sigma/md.cutoff)**6))), label="Shifted LJ potential")
        ax[0].plot(added, 0*added, color="orange")
        ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6))), label="LJ force")
        ax[1].plot(dist, ((4*md.epsilon/dist) * ((12*(md.sigma/dist)**12)-(6*(md.sigma/dist)**6)))-((4*md.epsilon/md.cutoff) * ((12*(md.sigma/md.cutoff)**12)-(6*(md.sigma/md.cutoff)**6))), label="Shifted LJ force")
        ax[1].plot(added, 0*added, color="orange")
    ax[0].legend()
    ax[1].set_ylim(top=100)
    ax[1].scatter(md.distancesContainer, md.forcesContainer, color="lime", s=1, label="All forces")
    ax[1].legend()
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.savefig(directory + "/potential.png", transparent=False, format="png")

def forces(md, directory: str) -> None:
    """
    Plot all the forces zoomed around cutoff.

    Parameters
    ----------
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    directory : str
        Directory in which to save the plot.

    Returns
    ----------
    A plot of all the zoomed forces around cutoff.
    """

    plt.clf()
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.axvline(x=md.cutoff, color="gray", linestyle="--")
    plt.scatter(md.distancesContainer, md.forcesContainer, color="orchid", s=4, label="Zoomed forces (all pairs)")
    plt.legend()
    plt.savefig(directory + "/continuous.png", transparent=False, format="png")

def sampleForces(md, directory: str) -> None:
    """
    Plot the forces of a sampled pair zoomed around cutoff.

    Parameters
    ----------
    md : MolecularDynamics
        An instance of the class MolecularDynamics.
    directory : str
        Directory in which to save the plot.

    Returns
    ----------
    A plot of the zoomed forces  of a sampled pair around cutoff.
    """

    plt.clf()
    plt.axhline(y=0, color="gray", linestyle="--")
    plt.plot(md.allforcesContainer, color="orchid", label="Zoomed forces over time (sampled pair)")
    plt.ylim(top=0.01, bottom=-0.001)
    plt.legend()
    plt.savefig(directory + "/continuous2.png", transparent=False, format="png")

#------------------------------------WORK-IN-PROGRESS---------------------------------------

#-----------------------------------------UNUSED--------------------------------------------

def densitySquares_cluster(directory, yDivision, perc, outputName, plot):

    with open(directory + os.sep +"classInstance.pkl", "rb") as f:
        md = pickle.load(f)

    step = int(np.shape(md.allPositions)[0]*perc/100)
    if step >= np.shape(md.allPositions)[0]: step = np.shape(md.allPositions)[0]-1
    positions = md.allPositions[step]

    sideSquares = md.box_size[1]/yDivision

    x = -md.box_size[0]/2
    y = -md.box_size[1]/2

    nParticles = np.zeros((yDivision * int(md.box_size[0] / md.box_size[1]), yDivision), dtype=int)

    for ii in range(yDivision):
        x = -md.box_size[0]/2
        for jj in range(yDivision * int(md.box_size[0] / md.box_size[1])):
            for particle in range(md.num_particles):
                if (x <= positions[particle, 0] < x + sideSquares) and (y <= positions[particle, 1] < y + sideSquares):
                    nParticles[jj, ii] += 1
            x += sideSquares
        y += sideSquares

    threshold = (np.max(nParticles) + np.min(nParticles))/2
    boxCheck = np.zeros((yDivision * int(md.box_size[0] / md.box_size[1]), yDivision), dtype=bool)

    y = -md.box_size[1]/2
    for ii in range(yDivision):
        x = -md.box_size[0]/2
        for jj in range(yDivision * int(md.box_size[0] / md.box_size[1])):
            if nParticles[jj, ii] > threshold: boxCheck[jj, ii] = True
            else: boxCheck[jj, ii] = False
            x += sideSquares
        y += sideSquares

    if plot:
        plt.figure()
        plt.imshow(np.fliplr(boxCheck.T), origin='upper', interpolation='nearest', extent=[-md.box_size[0]/2, md.box_size[0]/2, -md.box_size[1]/2, md.box_size[1]/2])
        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")
        plt.tight_layout()
        #plt.legend()
        plt.savefig(directory + f"/{outputName}.png", transparent=False, format="png")
    else: 
        # Complex network approach
        Nx, Ny = boxCheck.shape
        graph = nx.Graph()
        for ii in range(Nx):
            for jj in range(Ny):
                if boxCheck[ii,jj]:
                    graph.add_node((ii,jj))
                    for di, dj in [(-1,0),(1,0),(0,-1),(0,1)]:
                        ni = (ii + di) % Nx
                        nj = (jj + dj) % Ny
                        if boxCheck[ni,nj]:
                            graph.add_edge((ii,jj),(ni,nj))

        # Get all connected components
        components = list(nx.connected_components(graph))

        # Size of largest component
        giant_size = max(len(c) for c in components) if components else 0
        return giant_size/boxCheck.size
    
def gcc_time(directory, yDivision, outputName):

    percs = np.linspace(0, 100, 20)
    gccsize = np.zeros(np.shape(percs)[0])

    for ii, perc in enumerate(percs):
        gccsize[ii] = densitySquares_cluster(directory, yDivision, perc, outputName, plot=False)

    plt.figure()
    plt.plot(percs, gccsize, color="darkslategrey")
    plt.xlabel(r"perc")
    plt.ylabel(r"size GCC")
    plt.tight_layout()
    #plt.legend()
    plt.savefig(directory + f"/{outputName}.png", transparent=False, format="png")
