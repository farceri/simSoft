"""
mdFunctions 
----------
These are the functions that are needed to study the MolecularDynamics class.
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

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

#--------------------------------------MSD-SSF-ISF------------------------------------------

def compute_msd(positions):
    """Compute the Mean Squared Displacement."""
    initialPosition = positions[0]
    msd = np.zeros(np.shape(positions)[0]-1)
    for timestep in range(np.shape(positions)[0]-1):
        displacement = positions[timestep+1] - initialPosition
        msd[timestep] = np.mean(np.sum(displacement ** 2, axis=1))
    return msd

def compute_ssf(unwrappedPositions, num_particles, inf_lim, sup_lim, kMods):
    "Compute the Intermediate Scattering Function for different k modulus with mean over different t0."

    angles = np.arange(0, 2*np.pi, np.pi/4)
    ssf_int = np.zeros((np.shape(kMods)[0]), dtype=complex)
    mask = ~np.eye(num_particles, dtype=bool)
    ssf_total_self = np.ones((np.shape(kMods)[0]), dtype=complex)
    ssf_total_int = np.zeros((np.shape(kMods)[0]), dtype=complex)

    # Finding the maximum with the static structure factor
    for t0 in range(inf_lim, sup_lim): 
        ssf_int = np.zeros_like(ssf_int) 
        for ii, kMod in enumerate(kMods): 
            for angle in angles: 
                kVec = np.array((kMod * np.cos(angle), kMod * np.sin(angle))) 
                delta = unwrappedPositions[t0][:, None, :] - unwrappedPositions[t0] [None, :, :] 
                #delta -= self.box_size[0] * np.round(delta / self.box_size[0])
                delta = delta[mask].reshape(num_particles, num_particles-1, -1) 
                ssf_int[ii] += np.sum(np.exp(1j * np.tensordot(delta, kVec, axes=([2],[0]))))/(num_particles) 
            ssf_total_int[ii] += (ssf_int[ii]) / (np.shape(angles)[0])
        
    ssf_total_int /= (sup_lim - inf_lim)
    ssf_total_self = np.real(ssf_total_self)
    ssf_total_int = np.real(ssf_total_int)

    chosenk = np.array([kMods[np.argmax(ssf_total_self[:] + ssf_total_int[:])]])
    chosenk= chosenk[0]

    return ssf_total_self, ssf_total_int, chosenk

def compute_isf(unwrappedPositions, num_particles, inf_lim, sup_lim, chosenk):

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
        # Highlight one particle in red
        if ii != 0: scat_dic["scat{0}".format(ii)] = plt.scatter(positions[ii, 0, 0], positions[ii, 1, 0], 
                                                                 s=size,facecolor='cadetblue', edgecolor="black", linewidth=0.5)
        else: scat_dic["scat{0}".format(ii)] = plt.scatter(positions[ii, 0, 0], positions[ii, 1, 0], 
                                                           s=size, facecolor='firebrick', edgecolor="black", linewidth=0.5)
        #line_dic["line{0}".format(ii)] = plt.plot(positions[ii, 0, 0], positions[ii, 1, 0])[0]   

    # Gif visual setup    
    plt.xlim([-md.box_size[0]/2, md.box_size[0]/2])
    plt.ylim([-md.box_size[1]/2, md.box_size[1]/2])
    if md.interaction : plt.title(r"N=%d, T=%.1f, $\rho$=%.2f, $\sigma$=%.1f" %(md.num_particles, md.temperature, (md.num_particles*np.pi*((md.sigma/2)**2))/(md.box_size[0]**2), md.sigma))
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
