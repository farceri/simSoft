import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from numba import njit
import matplotlib.animation as animation
import sys

'''
Script to compute trajectories of random walkers and make histograms of the displacement pdf
Example usage: python3 randomWalk1D.py animate 10000 /path/to/figureName
'''

def makeAnimation(num_frames, traj, x, bins, mu, tau, dt, sigma, figureName):
    # This function makes an animation of the position pdf
    # Arguments: (num_frames) Number of frames in the animation
    # (traj) array of trajectories, typically numpy.array()
    # (x) plotting range, (bins) support of the histogram
    # (mu) mean of the distribution, (tau) correlation time
    # (dt) timestep, (sigma) noise amplitude

    # Set figure settings
    figs, axs = plt.subplots(figsize=(5,4), dpi=120)

    def update(frame):
        axs.clear()
        axs.set_xlim(bins[0], bins[-1])
        axs.set_ylim(0, 2.5)
        hist, _ = np.histogram(traj[:,frame + 1], bins=bins, density=True)
        axs.plot((bins[1:] + bins[:-1]) / 2, hist)
        y = norm.pdf(x, mu * (1 - np.exp(-(frame + 1) * dt / tau)), sigma**2 * np.sqrt(2) * (1 - np.exp(-2 * (frame + 1) * dt / tau)))
        axs.plot(x, y) 
        axs.set_title(f'Time: {(frame + 1) * dt:.4f}', fontsize=12)
        axs.set_xlim(bins[0], bins[-1])
        axs.set_ylim(0, 2.5)
        axs.set_xlabel("$Position,$ $x$", fontsize=14)
        axs.set_ylabel("$PDF(x)$", fontsize=14)
        plt.tight_layout()
        return axs.artists
    
    anim = animation.FuncAnimation(figs, update, frames=num_frames-1, interval=100, blit=False)
    if(figureName != 0):
        anim.save(f'{figureName}.mov', writer='ffmpeg', dpi=fig.dpi)
        print("Saving animation as", figureName)
    else:
        plt.show()


def plotSinglePDF(pos, x, bins, mu, tau, time, sigma):
    # This function plot the position pdf given a number of trials
    # Arguments: (x) plotting range, (bins) support of the histogram
    # (mu) mean of the distribution, (tau) correlation time
    # (time) current timestep, (sigma) noise amplitude

    # Compute pdf and plot it
    hist, edges = np.histogram(pos, bins=bins, density=True)
    centers = (edges[1:] + edges[:-1]) / 2
    ax.plot(centers, hist)
    y = norm.pdf(x, mu * (1 - np.exp(-time / tau)), sigma**2 * np.sqrt(2) * (1 - np.exp(-2 * time / tau)))
    ax.plot(x, y, color='r')


if __name__ == "__main__":
    # Global variables (sigma, mu, tau and dt can be inputs using sys.argv[])
    sigma = 1 / np.sqrt(2) # noise amplitude
    mu = 10 # mean of the distribution
    tau = 0.1 # correlation time
    dt = 1e-03 # time step
    total_time = 0.6 # total time
    num_steps = int(total_time / dt)
    print("Number of iterations:", num_steps)
    animate = sys.argv[1] # Flag for making animation
    ntrials = int(sys.argv[2]) # Number of realizations
    figureName = sys.argv[3]

    pos = np.zeros(ntrials)
    #MDS = np.zeros(ntrials)

    traj = np.zeros((ntrials,num_steps))
    #saveMDS = np.zeros((ntrials, num_steps))

    # Define bins for the histograms
    bins = np.linspace(-2, 14, 500)
    x = np.linspace(bins[0], bins[-1], 1000)
    fig, ax = plt.subplots(figsize=(5, 4), dpi=120)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 2.5)
    for i in range(num_steps):
        pos += (mu - pos) * dt / tau + sigma * np.sqrt(2. * dt / tau) * np.random.randn(ntrials)
        traj[:,i] = pos
        #MDS += (pos - mu)**2
        #saveMDS[:,i] = MDS
        if i in (5, 30, 90, 500):
            plotSinglePDF(pos, x, bins, mu, tau, i * dt, sigma)
    
    ax.set_xlabel("$Position,$ $x$", fontsize=14)
    ax.set_ylabel("$PDF(x)$", fontsize=14)
    plt.tight_layout()
    plt.pause(0.5)

    # Make the animation
    if animate == 'animate':
        print("Making animation")
        makeAnimation(num_steps, traj, x, bins, mu, tau, dt, sigma, figureName)
    else:
        plt.show()