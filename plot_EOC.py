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

from matplotlib.patches import Rectangle
# import matplotlib.patches as patches
import subprocess
# - - - - - - - - - - - - - - - - - -

figsize=15
fontsize=2*figsize
colors=['r','b']

# names of the files where you saved the output (EOC, etc)
filenames=['EOC_tmax10.dat','EOC_tmax100.dat']

fig, fig1 = plt.subplots(1,sharex=True,figsize=(figsize,figsize))

# read each file in filenames
for filename in filenames :
    # read the number of lines calling the terminal bash 'wc' command
    proc = subprocess.Popen(['wc',filename], stdout=subprocess.PIPE)
    tmp = proc.communicate()
    lines = int(tmp[0].split()[0]) - 1 # -1 to skip header
    N  =np.zeros(lines) # initiate arrays with the right size
    EOC=np.zeros(lines)
    f=open(filename) # open file
    tmp = f.readline() # read first line beforehand, to skip header
    for i in range(lines): # read actual lines with data
        tmp = f.readline() # read line
        N[i] = int(tmp.split()[1]) # store the integer corresponding to the 2nd element of this line
        EOC[i] = float(tmp.split()[4]) # store the float corresponding to the 5th element of this line
    fig1.plot(N,EOC,marker='o',markersize=fontsize,linewidth=0,color=colors[filenames.index(filename)]) # plot EOC function of N w/ round markers, no line and colors[index]
    # filenames.index(filename) gives the index of filename in the filenames list (0 for the 1st, 1 for the 2nd, etc)

expected_order_of_accuracy=1.
# plot a line @ expected_order_of_accuracy to indicate the limit the EOC markers should tend towards
fig1.plot(np.linspace(float(N[0])/2.,2*N[lines-1],2*lines),expected_order_of_accuracy*np.ones(2*lines),linewidth=2,color='k')
fig1.grid(linestyle='--',alpha=0.5)
fig1.minorticks_on()
fig1.set_xscale('log',basex=2)
fig1.set_yscale('log',basey=10)
fig1.set_xlim(15.,700.)
fig1.set_ylim(0.001,1.5)
fig1.tick_params(labelsize=fontsize)

# define legend and plot it in lower right corner
tmax=['tmax=10','tmax=100']
leg=2*['']
leg[0] = plt.Line2D([0],[0],marker='o',linewidth=0,color=colors[0],markerfacecolor=colors[0],markersize=fontsize/2)
leg[1] = plt.Line2D([0],[0],marker='o',linewidth=0,color=colors[1],markerfacecolor=colors[1],markersize=fontsize/2)
fig1.legend(leg,tmax,
	loc='lower right',shadow=True,fontsize=fontsize)

fig1.set_ylabel(r'EOC', fontweight='bold', fontsize=fontsize)
fig1.set_xlabel(r'Number of cells (log$_2$)',fontweight='bold',fontsize=fontsize)
fig.tight_layout()
plt.savefig('EOC.pdf',format='pdf',dpi=700,bbox_inches='tight')
plt.show()
plt.close()
