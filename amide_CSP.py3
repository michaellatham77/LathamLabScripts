#!/usr/bin/env python
#
# Script to make list for plotting amide chemical shift perturbations
# onto structure using pymol
#
# Based off of amide_CSP.py script
#
# MP Latham, Jan 17, 2018


import sys, re
import argparse
#import scipy as sp #7/16/21 SA
import numpy as np #7/16/21 SA
import matplotlib.pyplot as plt

#Argument parser
parser = argparse.ArgumentParser(
    prog = 'amide_CSP.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''\
    Script to calculated combined chemical shift perturbations.
    CSPs are weighted.
    
    Peaklist must be in the following format
    Assignment  w1  w2
    A30N-A30HN 120.2 8.23
    ''')
parser.add_argument('-apo', type=str, required=True, nargs=1, help='Sparky format amide peaklist of protein alone')
parser.add_argument('-bound', type=str, required=True, nargs=1, help='Sparky format amide peaklist of protein+ligand')
parser.add_argument('-pymol_bfactor', nargs=1, type=str, help='Make file for pymol data2bfactor and color_b. Include filename as argument')
parser.add_argument('--plot', action='store_true', help='Optional flag to plot bar chart of CSPs')

#Print help if no arguments are given
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

#Functions
def pymol_bfactor(amide_list, output):
    '''Make list for data2bfactor and color_b pymol scripts'''

    #Regular expression for 101IleH(orC)d1
    #p = re.compile("([0-9]*)([A-Z][a-z]*)([A-Za-z1-2]+)")

    #Regular expression for S17N-S17HN
    p = re.compile("([A-Z])([0-9]*)N-([A-Z])([0-9]*)HN")

    with open(output, 'w') as fo:

        for amide in amide_list:
            if '{' in amide[0] or '[' in amide[0]:
                pass
            else:
                m = p.findall(amide[0])

                if len(m) == 0:
                    '''This is to get around side chain assignments'''
                    pass
                else:
                    resi = m[0][1] #residue number
                    resn = m[0][0].upper() #residue type
                    #name = m[0][2].upper() #atom name
                    name = 'CA' #atom name

                    fo.write('{:4s} {:3s} {:3s} {:.4f}\n'.format(resi, resn, name, amide[3]))
                    print('{:4s} {:3s} {:3s} {:.4f}'.format(resi, resn, name, amide[3]))

def plot_csps(amide_list):
    '''make a bar chart of CSP vs residue number'''

    array = list()

    #Regular expression for S17N-S17HN
    p = re.compile("([A-Z])([0-9]*)N-([A-Z])([0-9]*)HN")

    for amide in amide_list:
        if '{' in amide[0]:
            pass
        else:
            m = p.findall(amide[0])

            if len(m) == 0:
                '''This is to get around side chain assignments'''
                pass
            else:
                array.append([int(m[0][1]), amide[3]])

    resi = [row[0] for row in array]
    csp = [row[1] for row in array]
    
    fig, ax = plt.subplots()
    graph = ax.bar(resi, csp)

    ax.set_xlabel('Residue Number')
    ax.set_ylabel('Chemical Shift Perturbation')

    plt.show()

#Read in data and put into Scipy record array
apo = np.genfromtxt(args.apo[0], comments='#', skip_header=1, dtype=None, usecols=(0,1,2),
			names=['Assign','15N','1H'])
bound = np.genfromtxt(args.bound[0], comments='#', skip_header=1, dtype=None, usecols=(0,1,2),
			names=['Assign','15N','1H'])

#Master list of chemical shift perturbations
csp_list = list()

for peak in apo:
    amide_group = peak['Assign']
    apo_H, apo_N = peak['1H'], peak['15N']

    if bound[np.where(bound['Assign'] == amide_group)]: #Find matches in the two lists
        bound_H = bound[np.where(bound['Assign'] == amide_group)]['1H'][0]
        bound_N = bound[np.where(bound['Assign'] == amide_group)]['15N'][0]

        delta_N2 = (apo_N - bound_N)**2.0
        delta_H2 = (apo_H - bound_H)**2.0
        csp = np.sqrt((delta_H2 + delta_N2 / 25.0) / 2.0)
        #csp = np.sqrt((apo_N - bound_N)**2.0 / 25.0 + (apo_H - bound_H)**2.0 / 2.0)
    	''' Garrett, DS et al (1997) Biochemistry, 36, 4393-8.'''

        csp_list.append((amide_group, (apo_H, apo_N), (bound_H, bound_N), csp))

    else:
        pass

if args.pymol_bfactor:
    pymol_bfactor(csp_list, args.pymol_bfactor[0])

if args.plot:
    plot_csps(csp_list)
