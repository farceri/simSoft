'''
Created by Francesco
March 31 2025
'''

#functions and scripts to visualize samples
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib import animation
from matplotlib import cm
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools
import sys
import os

def setAxes3D(ax):
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

def set3DPackingAxes(boxSize, ax):
    xBounds = np.array([0, boxSize[0]])
    yBounds = np.array([0, boxSize[1]])
    zBounds = np.array([0, boxSize[2]])
    ax.set_xlim(xBounds[0], xBounds[1])
    ax.set_ylim(yBounds[0], yBounds[1])
    ax.set_ylim(zBounds[0], zBounds[1])
    #ax.set_box_aspect(aspect = (1,1,1))
    #ax.set_aspect('equal', adjustable='box')
    setAxes3D(ax)

def plot3DPacking(dirName, figureName):
    boxSize = np.loadtxt(dirName + os.sep + 'boxSize.dat')
    rad = np.array(np.loadtxt(dirName + os.sep + 'rad.dat'))
    pos = np.array(np.loadtxt(dirName + os.sep + 'pos.dat'))
    pos[:,0] -= np.floor(pos[:,0]/boxSize[0]) * boxSize[0]
    pos[:,1] -= np.floor(pos[:,1]/boxSize[1]) * boxSize[1]
    pos[:,2] -= np.floor(pos[:,2]/boxSize[2]) * boxSize[2]
    fig = plt.figure(dpi=200)
    ax = Axes3D(fig)
    set3DPackingAxes(boxSize, ax)
    u = np.linspace(0, 2*np.pi, 120)
    v = np.linspace(0, np.pi, 120)
    colorId = getRadColorList(rad)
    for i in range(pos.shape[0]):
        x = pos[i,0] + rad[i]*np.outer(np.cos(u), np.sin(v))
        y = pos[i,1] + rad[i]*np.outer(np.sin(u), np.sin(v))
        z = pos[i,2] + rad[i]*np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x,y,z, color=colorId[i], rstride=4, cstride=4, alpha=1)
    plt.savefig('/home/francesco/Pictures/soft/packings/3d-' + figureName + '.png', transparent=True, format = 'png')
    plt.show()

def setAxes2D(ax):
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])

def setPackingAxes(boxSize, ax):
    if(boxSize.shape[0] == 1):
        ax.set_xlim(-boxSize, boxSize)
        ax.set_ylim(-boxSize, boxSize)
    else:
        ax.set_xlim(0, boxSize[0])
        ax.set_ylim(0, boxSize[1])
    ax.set_aspect('equal', adjustable='box')
    setAxes2D(ax)

def setBigBoxAxes(boxSize, ax, delta=1.1):
    if(boxSize.shape[0] == 1):
        ax.set_xlim(-boxSize, boxSize)
        ax.set_ylim(-boxSize, boxSize)
        xBounds = np.array([-boxSize*delta, boxSize*delta])
        yBounds = np.array([-boxSize*delta, boxSize*delta])
    else:
        xBounds = np.array([boxSize[0]*(1-delta), boxSize[0]*delta])
        yBounds = np.array([boxSize[1]*(1-delta), boxSize[1]*delta])
    ax.set_xlim(xBounds[0], xBounds[1])
    ax.set_ylim(yBounds[0], yBounds[1])
    ax.set_aspect('equal', adjustable='box')
    setAxes2D(ax)

def setInvisiblePackingAxes(boxSize, ax):
    setBigBoxAxes(boxSize, ax, 1.05)
    for spine in ax.spines.values():
        spine.set_visible(False)

def setCenteredPackingAxes(boxSize, frameSize, ax):
    # Centering the frame
    x_center = boxSize[0] / 2
    y_center = boxSize[1] / 2
    ax.set_xlim([x_center - frameSize[0] / 2, x_center + frameSize[0] / 2])
    ax.set_ylim([y_center - frameSize[1] / 2, y_center + frameSize[1] / 2])
    ax.set_aspect('equal', anchor='C')
    setAxes2D(ax)

def getRadColorList(rad):
    colorList = cm.get_cmap('viridis', rad.shape[0])
    colorId = np.zeros((rad.shape[0], 4))
    count = 0
    for particleId in np.argsort(rad):
        colorId[particleId] = colorList(count/rad.shape[0])
        count += 1
    return colorId

def getDoubleColorList(rad, numA=0):
    colorId = np.zeros((rad.shape[0], 4))
    colorId[:numA] = [0,1,0,0.2]
    colorId[numA:] = [0,0,1,0.2]
    return colorId

def getPBCPositions(fileName, boxSize):
    pos = np.array(np.loadtxt(fileName), dtype=np.float64)
    pos[:,0] -= np.floor(pos[:,0]/boxSize[0]) * boxSize[0]
    pos[:,1] -= np.floor(pos[:,1]/boxSize[1]) * boxSize[1]
    return pos

def shiftPositions(pos, boxSize, xshift, yshift):
    pos[:,0] += xshift
    pos[:,1] += yshift
    pos[:,0] -= np.floor(pos[:,0]/boxSize[0]) * boxSize[0]
    pos[:,1] -= np.floor(pos[:,1]/boxSize[1]) * boxSize[1]
    return pos

def getDirSep(dirName, fileName):
    if(os.path.exists(dirName + os.sep + fileName + ".dat")):
        return "/"
    else:
        return "/../"

def plotPacking(dirName, figureName, quiver=False, lj=False, shiftx=0, shifty=0, double=False, numA=0, alpha=0.6, lw=0.3):
    sep = getDirSep(dirName, 'boxSize')
    boxSize = np.atleast_1d(np.loadtxt(dirName + sep + 'boxSize.dat'))
    pos = getPBCPositions(dirName + os.sep + 'pos.dat', boxSize)
    rad = np.array(np.loadtxt(dirName + sep + 'rad.dat'))
    # plot particle size as interaction size
    if lj: rad *= 2**(1/6)
    # apply shift if different from zero
    if(shiftx != 0 or shifty != 0):
        pos = shiftPositions(pos, boxSize, shiftx, shifty)
    # make figure
    fig, ax = plt.subplots(dpi=200)
    setPackingAxes(boxSize, ax)
    if double:
        colorId = getDoubleColorList(rad, numA)
    else:
        colorId = getRadColorList(rad)
    if quiver:
        vel = np.array(np.loadtxt(dirName + os.sep + 'vel.dat'))
    for particleId in range(rad.shape[0]):
        x = pos[particleId,0]
        y = pos[particleId,1]
        r = rad[particleId]
        #print('particle', particleId, 'position:', x, y)
        if quiver:
            ax.add_artist(plt.Circle([x, y], r, edgecolor='k', facecolor=colorId[particleId], alpha=alpha, linewidth=lw))
            vx = vel[particleId,0]
            vy = vel[particleId,1]
            ax.quiver(x, y, vx, vy, facecolor='k', linewidth=0.1, width=0.001, scale=80, headlength=5, headaxislength=5, headwidth=5, alpha=0.6)
        else:
            ax.add_artist(plt.Circle([x, y], r, edgecolor='k', facecolor=colorId[particleId], alpha=alpha, linewidth=lw))    
        #label = ax.annotate(str(particleId), xy=(x, y), fontsize=4, verticalalignment='center', horizontalalignment='center')
    if quiver:
        figureName = '/home/francesco/Pictures/soft/packings/velmap-' + figureName + '.png'
    else:
        figureName = '/home/francesco/Pictures/soft/packings/' + figureName + '.png'
    plt.tight_layout()
    plt.savefig(figureName, transparent=True, format='png')
    plt.show()

########################## check energy conservation ##########################
def plotEnergy(dirName, figureName):
    if(os.path.exists(dirName + os.sep + "energy.dat")):
        energy = np.loadtxt(dirName + os.sep + "energy.dat")
        print("potential energy:", np.mean(energy[:,2]), "+-", np.std(energy[:,2]), "  std/mean:", np.std(energy[:,2])/np.abs(np.mean(energy[:,2])))
        print("temperature:", np.mean(energy[:,3]), "+-", np.std(energy[:,3]), "  std/mean:", np.std(energy[:,3])/np.abs(np.mean(energy[:,3])))
        fig, ax = plt.subplots(figsize=(6,4.5), dpi = 120)
        ax.plot(energy[:,0], energy[:,2], linewidth=1.5, color='k', linestyle='solid', label="$U$")
        ax.plot(energy[:,0], energy[:,3], linewidth=1.5, color='r', linestyle='dashed', label="$K$")
        ax.plot(energy[:,0], energy[:,4], linewidth=4, color='b', linestyle='solid', alpha=0.3)
        ax.plot(energy[:,0], energy[:,4], linewidth=1.5, color='b', linestyle='dotted', label="$E$")
        ax.set_ylabel("$Energy$", fontsize=16)
        print("\ntotal energy per particle:", np.mean(energy[:,4]), "+-", np.std(energy[:,4]), "  std/mean:", np.std(energy[:,4])/np.abs(np.mean(energy[:,4])))
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xlabel("$Simulation$ $step$", fontsize=16)
        ax.legend(fontsize=12, loc='upper left')
        plt.tight_layout()
        figureName = "/home/francesco/Pictures/soft/energy-" + figureName
        fig.savefig(figureName + ".png", transparent=True, format = "png")
        plt.show()
    else:
        print("no energy.dat file was found in", dirName)

def readFromParams(dirName, paramName):
    name = None
    with open(dirName + os.sep + "params.dat") as file:
        for line in file:
            name, scalarString = line.strip().split("\t")
            if(name == paramName):
                return float(scalarString)
    if(name == None):
        print("The variable", paramName, "is not saved in this file")
        return None
    
def getOrderedDirectories(dirName):
    listDir = []
    listScalar = []
    for dir in os.listdir(dirName):
        if(os.path.isdir(dirName + os.sep + dir)):
            listDir.append(dir)
            listScalar.append(dir.strip('t'))
    listScalar = np.array(listScalar, dtype=np.int64)
    listDir = np.array(listDir)
    listDir = listDir[np.argsort(listScalar)]
    listScalar = np.sort(listScalar)
    return listDir, listScalar

def plotTest2Forces(dirName, which=None):
    timeStep = readFromParams(dirName, "dt")
    dirList, timeList = getOrderedDirectories(dirName)
    timeList = np.array(timeList) * timeStep
    force0x = []
    force0y = []
    force1x = []
    force1y = []
    for d in range(dirList.shape[0]):
        force = np.loadtxt(dirName + os.sep + dirList[d] + "/forces.dat")
        force0x.append(force[0,0])
        force0y.append(force[0,1])
        force1x.append(force[1,0])
        force1y.append(force[1,1])
    fig = plt.figure(figsize = (7, 5), dpi = 120)
    ax = fig.gca()
    if(which=='0'):
        ax.plot(timeList, force0x, linewidth=1, color='k', marker='o', fillstyle='none', label="$p0, x$")
        ax.plot(timeList, force0y, linewidth=1, color='k', marker='v', fillstyle='none', label="$p0, y$")
    elif(which=='1'):
        ax.plot(timeList, force1x, linewidth=1, color='g', marker='o', fillstyle='none', label="$p1, x$")
        ax.plot(timeList, force1y, linewidth=1, color='g', marker='v', fillstyle='none', label="$p1, y$")
    else:
        ax.plot(timeList, force0x, linewidth=1, color='k', marker='o', fillstyle='none', label="$p0, x$")
        ax.plot(timeList, force0y, linewidth=1, color='k', marker='v', fillstyle='none', label="$p0, y$")
        ax.plot(timeList, force1x, linewidth=1, color='g', marker='o', fillstyle='none', label="$p1, x$")
        ax.plot(timeList, force1y, linewidth=1, color='g', marker='v', fillstyle='none', label="$p1, y$")
    ax.legend(fontsize=10, loc='best')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel("$Simulation$ $step$", fontsize=15)
    ax.set_ylabel("$Forces$", fontsize=15)
    plt.tight_layout()
    plt.show()

# Main script with plotting options
if __name__ == '__main__':
    dirName = sys.argv[1]
    whichPlot = sys.argv[2]
    figureName = sys.argv[3]

    if(whichPlot == 'plot'):
        plotPacking(dirName, figureName, shiftx=float(sys.argv[4]), shifty=float(sys.argv[5]))
    
    if(whichPlot == 'plotvel'):
        plotPacking(dirName, figureName, quiver=True, shiftx=float(sys.argv[4]), shifty=float(sys.argv[5]))

    elif(whichPlot == 'plotlj'):
        plotPacking(dirName, figureName, lj=True, shiftx=float(sys.argv[4]), shifty=float(sys.argv[5]))

    elif(whichPlot == 'plotljvel'):
        plotPacking(dirName, figureName, quiver=True, lj=True, shiftx=float(sys.argv[4]), shifty=float(sys.argv[5]))

    elif(whichPlot == 'plot2'):
        plotPacking(dirName, figureName, double=True, numA=int(sys.argv[4]), shiftx=float(sys.argv[5]), shifty=float(sys.argv[6]))

    elif(whichPlot == 'plot2vel'):
        plotPacking(dirName, figureName, double=True, numA=int(sys.argv[4]), quiver=True, shiftx=float(sys.argv[5]), shifty=float(sys.argv[6]))

    elif(whichPlot == 'plot2lj'):
        plotPacking(dirName, figureName, double=True, lj=True, numA=int(sys.argv[4]), shiftx=float(sys.argv[5]), shifty=float(sys.argv[6]))

    elif(whichPlot == 'plot2ljvel'):
        plotPacking(dirName, figureName, double=True, lj=True, numA=int(sys.argv[4]), quiver=True, shiftx=float(sys.argv[5]), shifty=float(sys.argv[6]))

    elif(whichPlot == 'energy'):
        plotEnergy(dirName, figureName)

    elif(whichPlot == 'test2'):
        plotTest2Forces(dirName, which=sys.argv[4])

    else:
        print('Please specify the type of plot you want')