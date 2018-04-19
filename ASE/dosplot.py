#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from pylab import *

import os, sys, pickle

#Read input pickle files
def read_dos(dir,tetrahedra=False):
    import pickle
    try:
        if tetrahedra==True:
            f = open(dir + '/dos_tetra.pickle')
            energies, dos, pdos = pickle.load(f)
            f.close()
        else:
            f = open(dir + '/dos.pickle')
            energies, dos, pdos = pickle.load(f)
            f.close()
    except:
        print "No Density of States DATA Found."
        sys.exit(1)
    return energies, dos, pdos

rcParams['figure.figsize'] = 6*1.67323,4*1.67323
rcParams['ps.useafm'] = True
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

rcParams['pdf.fonttype'] = 42
matplotlib.rc('xtick.major', size=6)
matplotlib.rc('xtick.minor', size=3)
matplotlib.rc('ytick.major', size=6)
matplotlib.rc('ytick.minor', size=3)
matplotlib.rc('lines', markeredgewidth=0.5*2)
matplotlib.rc('font', size=12*2)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('PDOS')

num1 = 00 # Atom index for which the DOS will be extracted and plotted

energies,dos,pdos = read_dos('./')
resov_dos1 = pdos[num1]['d'][0] #Plot the projected density of the 'd' states 

ax.plot(energies,resov_dos1,
        color='k',
        linestyle='-',
	label= 'M d DOS') #Metal d-states (legend)

ax.axis([-10.,10.,0,10.])
ax.minorticks_on()


leg=plt.legend(loc=1,prop={'size':24},numpoints=1)
leg.draw_frame(False)

left  = 0.125  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.2   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.35   # the amount of width reserved for blank space between subplots
hspace = 0.35   # the amount of height reserved for white space between subplots
subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

plt.savefig('xyz.png', format='png')
plt.show()
