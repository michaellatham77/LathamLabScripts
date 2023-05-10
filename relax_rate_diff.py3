#!/usr/bin/env python3
#
# Script to take the difference between two sets of relaxation rates
# Script will also plot the difference vs residue number and/or
# make a table for pymol
#
# Based off of amide_CSP.py script
#
# MP Latham, April 21, 2019


import sys, re
import argparse
import numpy as np
import matplotlib.pyplot as plt

#Argument parser
parser = argparse.ArgumentParser(
    prog = 'relax_rate_diff.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''\
    Script to take the difference of relaxation rates
    with optional plotting of difference vs residue number
    or table for pymol
    
    Assumes that tables of relaxation rates are from exp_fit_pipe.py or similar
    and have peak names in column 1, rates in column 4, and errors in column 5.
    # Peak Name      Amplitude  Esd(Amp.)       Rate  Esd(Rate)
    V129N-V129HN     115730.13     523.85   17.59608    0.20596
    D128N-D128HN     117140.80     710.17   19.85383    0.28881
    ''')
parser.add_argument('-rate1', type=str, required=True, nargs=1, help='Table of first relaxation rates')
parser.add_argument('-rate2', type=str, required=True, nargs=1, help='Table of second relaxation rates')
parser.add_argument('-pymol_bfactor', nargs=1, type=str, help='Make file for pymol data2bfactor and color_b. Include filename as argument')
parser.add_argument('--plot', action='store_true', help='Optional flag to plot delta(rate) vs residue number')

#Print help if no arguments are given
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

#Functions
def pymol_bfactor(amide_list, output):
    '''Make list for data2bfactor and color_b pymol scripts'''

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

                    fo.write('{:4s} {:3s} {:3s} {:.4f}\n'.format(resi, resn, name, amide[1]))
                    print('{:4s} {:3s} {:3s} {:.4f}'.format(resi, resn, name, amide[1]))

def plot_rate(amide_list):
    '''make a plot of delta rate vs residue number'''

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
                array.append([int(m[0][1]), amide[1], amide[2]])

    resi = [row[0] for row in array]
    rate = [row[1] for row in array]
    err = [row[2] for row in array]
    
    fig, ax = plt.subplots()
    graph = ax.errorbar(resi, rate, yerr=err, xerr=None,
                fmt='o', markeredgewidth=1.0,
                linewidth=2.0, barsabove=False,
                markersize=8)

    ax.set_xlabel('Residue Number')
    ax.set_ylabel('Delta(R)')

    plt.show()

#Read in data and put into numpy record array
rate1 = np.genfromtxt(args.rate1[0], comments='#', skip_header=1, dtype=None, usecols=(0,3,4),
			names=['Assign','rate','err'], encoding=None)
rate2 = np.genfromtxt(args.rate2[0], comments='#', skip_header=1, dtype=None, usecols=(0,3,4),
			names=['Assign','rate','err'], encoding=None)

#Master list of relaxation rate differences 
rate_list = list()

for peak in rate1:
    amide_group = peak['Assign']
    rate_1, err_1 = peak['rate'], peak['err']

    if rate2[np.where(rate2['Assign'] == amide_group)]: #Find matches in the two lists
        rate_2 = rate2[np.where(rate2['Assign'] == amide_group)]['rate'][0]
        err_2 = rate2[np.where(rate2['Assign'] == amide_group)]['err'][0]

        delta_rate = rate_1 - rate_2
        delta_err = np.abs(delta_rate * np.sqrt((err_1/rate_1)**2.0 + (err_2/rate_2)**2.0))

        rate_list.append((amide_group, delta_rate, delta_err))

        if not args.pymol_bfactor:
            print('{:s} {:.2f} {:.2f}'.format(amide_group, delta_rate, delta_err))

    else:
        pass

if args.pymol_bfactor:
    pymol_bfactor(rate_list, args.pymol_bfactor[0])

if args.plot:
    plot_rate(rate_list)
