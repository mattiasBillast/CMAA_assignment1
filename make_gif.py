import numpy as np
from numpy import linspace
import math
# import time as tm
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import csv
import scipy.stats
from scipy.stats import norm
#from scipy.stats import skewnorm
from scipy.interpolate import interp2d, interp1d
import scipy.ndimage
from scipy import signal
import heapq
# import matplotlib.patches as patches

import imageio

from matplotlib.patches import Rectangle
# import matplotlib.patches as patches
import subprocess
# - - - - - - - - - - - - - - - - - -

figsize=15
fontsize=2*figsize

# TO BE SET BY HAND
iniCond='linear_acoustic'
Nfiles=31

filename_base=iniCond+'_'
nw=3 # number of variables (mass density, 3 velocity components and pressure)

# where the snapshot are stored before being bound together into a GIF
images = []

# read each file in filenames
for i in range(Nfiles):
    print 'Processing file ', i, 'in', Nfiles
    filename=filename_base+"%04d" % (i,)+'.dat' # "%04d" % (i,) to reproduce front zeros padding
    fig, fig1 = plt.subplots(nrows=1,ncols=3,sharex=True,figsize=(3*figsize,figsize))
    # read the number of lines calling the terminal bash 'wc' command
    proc = subprocess.Popen(['wc',filename], stdout=subprocess.PIPE)
    tmp = proc.communicate()
    lines = int(tmp[0].split()[0]) # no header this time
    x=np.zeros(lines) # initiate arrays with the right size
    w=np.zeros((lines,nw)) # 2D array
    f=open(filename) # open file
    for k in range(lines): # read actual lines with data
        tmp = f.readline() # read line
        x[k] = float(tmp.split()[0]) # store the float corresponding to the 1st element of this line (position)
        for j in range(nw): # starts @ 0, stops @ nw-1 included
            w[k][j] = float(tmp.split()[1+j]) # store the float corresponding to the 5th element of this line
    # plot only mass density, vx and pressure
    fig1[0].plot(x,w[:,0],marker='o',markersize=fontsize/5,linewidth=0,color='k') # plot EOC function of N w/ round markers, no line and colors[index]
    fig1[1].plot(x,w[:,1],marker='o',markersize=fontsize/5,linewidth=0,color='k') # plot EOC function of N w/ round markers, no line and colors[index]
    fig1[2].plot(x,w[:,4],marker='o',markersize=fontsize/5,linewidth=0,color='k') # plot EOC function of N w/ round markers, no line and colors[index]
    # beware, manually set domain range...
    if (iniCond is 'shock_tube'):
        # for density
        fig1[0].set_ylim(1.,8.)
        # for speed
        fig1[1].set_ylim(0.,1.)
        # for pressure
        fig1[2].set_ylim(1.,8.)
    elif (iniCond is 'linear_acoustic'):
        # for density
        fig1[0].set_ylim(0.992,1.008)
        # fig1[0].set_ylim(0.6,1.4) # for NL, with 0.07 amplitude on 0.1 pressure level
        # for speed
        fig1[1].set_ylim(-0.0015,0.0015)
        # fig1[1].set_ylim(-0.1,0.1) # for NL
        # for pressure
        fig1[2].set_ylim(0.1,0.1011)
        # fig1[2].set_ylim(0.09,0.18) # for NL
    # write the index of the snapshot
    fig1[2].text(0.9, 0.9,str(i),horizontalalignment='center',verticalalignment='center',transform=fig1[2].transAxes,fontsize=2*fontsize)
    fig1[0].set_ylabel(r'$\rho$', fontweight='bold', fontsize=fontsize)
    fig1[0].set_xlabel(r'x',fontweight='bold',fontsize=fontsize)
    fig1[1].set_ylabel(r'v', fontweight='bold', fontsize=fontsize)
    fig1[1].set_xlabel(r'x',fontweight='bold',fontsize=fontsize)
    fig1[2].set_ylabel(r'P', fontweight='bold', fontsize=fontsize)
    fig1[2].set_xlabel(r'x',fontweight='bold',fontsize=fontsize)
    for ax in fig1.flat:
        ax.set_xlim(1.1*x[0],1.1*x[lines-1])
        ax.grid(linestyle='--',alpha=0.5)
        ax.minorticks_on()
        ax.tick_params(labelsize=fontsize)
    fig.tight_layout()
    plt.savefig(iniCond+'.png',format='png',dpi=70,bbox_inches='tight') # overwrite same file each time
    images.append(imageio.imread(iniCond+'.png')) # but before, add this file to images

imageio.mimsave(iniCond+'.gif',images,duration=30./float(Nfiles)) # duration per snapshot such as the whole GIF lasts 10s
