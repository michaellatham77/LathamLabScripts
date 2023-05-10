#!/usr/bin/env python3
#
# Script to plot PALES output Calculated vs Experimental RDCs
#
# MP Latham, July 11, 2017
#

import sys, re
import argparse

import nmrglue as ng
import matplotlib.pyplot as plt

from numpy import corrcoef, sqrt, mean, amin, amax
#from matplotlib.pyplot import connect

try:
    from tooltip_v3 import *
except ImportError:
    print('\ntooltip library is not in path.\n')
    pass

# Command line arguments
parser = argparse.ArgumentParser(description='Make PALES RDC input from difTab.tcl output and pdb file.')
parser.add_argument('-pales', type=str, nargs=1, help='pales output file')
parser.add_argument('-pdf', type=str, nargs=1, help='pdf output file') # SA 8/13/22

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

# Read in RDC table
pcomment, pformat, ptable = ng.fileio.pipe.read_table(args.pales[0])

pearsons = corrcoef(ptable['D_OBS'],ptable['D'])[0,1]
rmsd = sqrt(mean((ptable['D_OBS']-ptable['D'])**2.0))

fig = plt.figure(linewidth=1.0)
ax = fig.add_subplot(111)

#ax.scatter(ptable['D_OBS'],ptable['D'])
ax.errorbar(ptable['D_OBS'], ptable['D'], xerr=ptable['DD'], yerr=None,
            fmt='bo', linewidth=1.0, 
            markerfacecolor='w',
            markeredgewidth=1.0,
            markeredgecolor='b',
            barsabove=False)

lims = [amin([ax.get_xlim(), ax.get_ylim()]),
        amax([ax.get_xlim(), ax.get_ylim()]),]

ax.plot(lims, lims, 'k-', alpha=0.75)

ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)

ax.set_xlabel('Observed RDCs (Hz)')
ax.set_ylabel('Calculated RDCs (Hz)')

if pearsons > 0:
    xpos = 0.05
    ypos = 0.95
    ha= 'left'
else:
    xpos = 0.95
    ypos = 0.5
    ha = 'right'

ax.text(xpos, ypos, 'rmsd = {:0.3f} Hz\nR = {:0.3f}'.
        format(rmsd, pearsons),
        horizontalalignment=ha,
        verticalalignment='top',
        transform=ax.transAxes,
        fontsize=16)

try:
    af = AnnoteFinder(ptable['D_OBS'], ptable['D'], ptable['RESID_I'])
    plt.connect('button_press_event', af)
except NameError:
    pass

# SA 8/13/22
if args.pdf:
    plt.savefig(args.pdf[0])
plt.show()
