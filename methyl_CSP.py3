#!/usr/bin/env python3

import sys, re
import argparse
import numpy as sp

#Argument parser
parser = argparse.ArgumentParser(
    prog = 'methyl_CSP.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''\
    Script to calculated combined chemical shift perturbations.
    CSPs are weighted based on residue specific standard deviations of bmrb values.

    Peaklist must be in the following format
    Assignment  w1  w2
    70IleHd1 13.89 1.01
    ''')
parser.add_argument('-apo', type=str, required=True, nargs=1, help='Sparky format methyl peaklist of protein alone')
parser.add_argument('-bound', type=str, required=True, nargs=1, help='Sparky format methyl peaklist of protein+ligand')
parser.add_argument('-pymol_bfactor', action='store', type=str, help='Make file for pymol data2bfactor and color_b. Include filename as argument')

#Print help if no arguments are given
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

#Functions
def combine_cs(assign, h_cs, c_cs):
    '''Calculate methyl group weighted combined chemical shifts.
       Weights are standard deviate reported in BMRB.'''

    if 'ile' in assign.lower():
        h_weight = 0.29
        c_weight = 1.65
    elif 'leu' in assign.lower():
        h_weight = 0.28
        c_weight = 1.6
    elif 'val' in assign.lower():
        h_weight = 0.27
        c_weight = 1.4
    elif 'met' in assign.lower():
        h_weight = 0.41
        c_weight = 1.54
    else:
        #if unassigned use the average of ile, leu, and val
        h_weight = 0.28
        c_weight = 1.55

    return h_cs / h_weight, c_cs/ c_weight

def pymol_bfactor(methyl_list, output):
    '''Make list for data2bfactor and color_b pymol scripts'''

    csps = []

    #Regular expression for 101IleH(orC)d1
    p = re.compile("([0-9]*)([A-Z][a-z]*)([A-Za-z1-2]+)")

    #Regular expression for sparky format M1CE-M1QE
    q = re.compile("([A-Z])([0-9]*)([CA-Z1-2]+)")

    with open(output, 'w') as fo:

        for methyl in methyl_list:
            if '-' in methyl[0]:
                m = q.findall(methyl[0])

                resi = m[0][1] #residue number

                if m[0][0] == 'I':
                    resn = 'ILE'
                elif m[0][0] == 'L':
                    resn = 'LEU'
                elif m[0][0] == 'V':
                    resn = 'VAL'
                elif m[0][0] == 'M':
                    resn = 'MET'
                else:
                    pass
		
                name = m[0][2].upper()

            else:
                m = p.findall(methyl[0])

                resi = m[0][0] #residue number
                resn = m[0][1].upper() #residue type
                name = m[0][2].upper().replace('H','C') #atom name, replace an H with C

            fo.write('{:4s} {:3s} {:3s} {:.4f}\n'.format(resi, resn, name, methyl[3]))
            print('{:4s} {:3s} {:3s} {:.4f}'.format(resi, resn, name, methyl[3]))

            csps.append(methyl[3])

    print('\n')
    print('The average CSP is {:0.3f}'.format(sp.mean(csps)))
    print('The standard deviation of CSPs is {:0.3f}'.format(sp.std(csps, ddof=1)))

    return

#Read in data and put into Scipy record array
apo = sp.genfromtxt(args.apo[0], comments='#', skip_header=1, usecols=(0,1,2), dtype=None,
                    names=['Assign','13C','1H'], encoding=None)
bound = sp.genfromtxt(args.bound[0], comments='#', skip_header=1, usecols=(0,1,2), dtype=None,
                    names=['Assign','13C','1H'], encoding=None)

#Master list of chemical shift perturbations
csp_list = list()

for peak in apo:
    methyl_group = peak['Assign']

    if '?' in methyl_group or '[' in methyl_group or '{' in methyl_group or 'None' in methyl_group:
        pass

    else:
        apo_H, apo_C = combine_cs(methyl_group, peak['1H'], peak['13C'])

        #if bound[sp.where(bound['Assign'] == methyl_group)]: #Find matches in the two lists
        try:
            bound_1H = bound[sp.where(bound['Assign'] == methyl_group)]['1H'][0]
            bound_13C = bound[sp.where(bound['Assign'] == methyl_group)]['13C'][0]

            bound_H, bound_C = combine_cs(methyl_group, bound_1H, bound_13C)
            csp = sp.sqrt((apo_H - bound_H)**2.0 + (apo_C - bound_C)**2.0)

            csp_list.append((methyl_group, (apo_H, apo_C), (bound_H, bound_C), csp))

            if not args.pymol_bfactor:
                print('{:s} {:.3f}'.format(methyl_group, csp))

        except ValueError:
            print('There is something wrong with {:s}.'.format(bound['Assign']))

        except:
            pass

if args.pymol_bfactor:
    pymol_bfactor(csp_list, args.pymol_bfactor)
else:
    print('\nThe average CSP is {:.3f}.'.format(sp.mean([row[3] for row in csp_list])))
    print('The standard deviation of  CSPs is {:.3f}.'.format(sp.std([row[3] for row in csp_list], ddof=1)))
