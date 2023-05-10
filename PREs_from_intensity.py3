#!/usr/bin/env python3
#
# Script to calculate PREs from intensity ratios
# Intensities come from CCPN analysis peak lists
#
# Can normalize ratio to a given peak in the reference list
#
# Can make output tables for pymol and HADDOCK
# For HADDOCK, can include the two segids and an offset of residue numbers
# as additional inputs
#
# MP Latham, Sept. 23, 2022
#

import sys, re
import argparse

import numpy as np

from scipy.optimize import fsolve
from scipy.constants import mu_0, pi

parser = argparse.ArgumentParser(description='Script to take ratio of intensities from two peak lists')
parser.add_argument('-ref_list', type=str, required=True, help='peak list for reference spectrum')
parser.add_argument('-exp_list', type=str, required=True, help='peak list for experimental spectrum')
parser.add_argument('-ref_peak', type=str, help='peak to normalize intensities to')
parser.add_argument('-frac_label', type=float, default=1.0, help='fractional occupancy of sping label')
parser.add_argument('-pymol_bfactor', nargs=1, type=str,
                    help='make file for pymol data2bfactor and color_b. Include filename as argument')
parser.add_argument('-haddock', nargs=2,
                    help='calculate distances based on intensity ratio for HADDOCK. Include residue number where spin label is attached and filename as argument')
parser.add_argument('-segid', nargs=2, type=str,
                    help='two segids to pass to haddock. The first is the PRE residue, and the second is the molecule with spin label')
parser.add_argument('-offset', type=int, default=0, help='number to add to PRE residues')
parser.add_argument('-cutoff', type=float, default=0.33, help='Intensity ratio to consider for PRE distances')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

#Functions
def return_methyl_name(assign):
    '''Function to get methyl information from analysis assignment'''

    #Regular expression for 101IleH(orC)d1
    p = re.compile("([0-9]*)([A-Z][a-z]*)([A-Za-z1-2]+)") 
    #Regular expression for 101Ile[10]
    q = re.compile("([0-9]*)([A-Za-z]*)\[[^\]]*\]") 

    if '{' in assign or '?' in assign or 'None' in assign:
        resi = 'None'
        resn = 'None'
        name = 'None'

    elif '[' in assign:
        m = q.findall(assign)

        resi = m[0][0] #residue number
        resn = m[0][1].upper() #residue type
        name = 'None' #Don't know the atom name

    else:
        m = p.findall(assign)

        resi = m[0][0] #residue number
        resn = m[0][1].upper() #residue type
        name = m[0][2].upper().replace('H','C') #residue name

    return resi, resn, name

def R2_from_linewidth(width):
    '''Estimate R2 from linewidth'''

    return width * pi

def calc_pre_from_inten_ratio(r2, ratio):
    '''Calculate PRE from intensity ratio based on equation from Battiste and
       Wagner (2000) Biochemistry (eq 4).
       Time hard coded for HMQC taua = 0.00177 sec.
       Using fsolve to find the root of the equation.'''

    func = lambda pre: ratio - (r2 * np.exp(-pre * 4. * 0.00177)) / (r2 + pre)
    PRE = fsolve(func, 10)

    return PRE[0]

def calc_dist(pre, tauc, B0):
    '''Calculate distances from PRE'''

    gamma_H = 26.7522128e+07 #1H gyromagnetic ratio (T^-1 s^-1)
    g = -2.0023 #electron g value
    beta = 9.27401e-24 #Bohr magneton (J T^-1)
    #spin = 7./2. #Gd3+ spin
    spin = 1./2. #nitroxide spin label

    C = (mu_0**2) * (gamma_H**2) * (g**2) * (beta**2) * spin*(spin + 1)
    C = C/(15 * (4 * pi)**2)
    
    w_h = 2 * pi * B0
    B = (4.0*tauc + 3.0*tauc / (1 + (tauc * w_h)**2))

    dist = ((C / pre) * B)**(1./6.)

    #Return distance in Angstroms
    return dist*1e10


#Read in reference and experimental peak lists
ref_list = np.genfromtxt(args.ref_list, dtype=None, encoding=None, comments='#',
                         skip_header=1, names=['F1','F2','assign','height','width'])
exp_list = np.genfromtxt(args.exp_list, dtype=None, encoding=None, comments='#',
                         skip_header=1, names=['F1','F2','assign','height','width'])

#Deal with scaling factor
if args.ref_peak:
    if ref_list[np.where(ref_list['assign'] == args.ref_peak)].size > 0:
        ref_int = ref_list[np.where(ref_list['assign'] == args.ref_peak)]['height']
    else:
        print('Reference peak {:s} does not exist in reference list {:s}.\n'.
                format(args.ref_peak, args.ref_list))
        sys.exit()

    if exp_list[np.where(exp_list['assign'] == args.ref_peak)].size > 0:
        exp_int = exp_list[np.where(exp_list['assign'] == args.ref_peak)]['height']
    else:
        print('Reference peak {:s} does not exist in experiment list {:s}.\n'.
                format(args.ref_peak, args.exp_list))
        sys.exit()
    
    scale_factor = ref_int / exp_int
else:
    scale_factor = 1.0

#Loop through reference list
data_table = list()
for peak in ref_list:
    resi, resn, name = return_methyl_name(peak['assign'])
    if resi == 'None':
        pass
    else:
        ref_peak_int = peak['height']

        if exp_list[np.where(peak['assign'] == exp_list['assign'])].size > 0:
            exp_peak_int = exp_list[np.where(peak['assign'] == exp_list['assign'])]['height'] 

            #Deal with fractional labeling here
            int_ratio = scale_factor * (exp_peak_int / ref_peak_int)
            data_table.append([peak['assign'], resi, resn, name, int_ratio[0]])

#Make outputs
if args.pymol_bfactor:
    '''Make text file for data2bfactor and color_b pymol scripts'''
    
    if args.offset:
        offset = args.offset[0]
    else:
        offset = 0

    with open(args.pymol_bfactor[0], 'w') as fo:
        for residue in data_table:
            resi, resn, name = residue[1], residue[2], residue[3]
            resi = str(offset + int(resi))
            int_ratio = residue[4]

            fo.write('{:4s} {:3s} {:3s} {:.3f}\n'.format(resi, resn, name, int_ratio))
            print('{:4s} {:3s} {:3s} {:.3f}'.format(resi, resn, name, int_ratio))

if args.haddock:
    '''Make text file for HADDOCK unambiguous restraints'''
    #Sort output file for haddock and the residue number where the spin label was attached
    try:
        int(args.haddock[0])
    except:
        haddock_out = args.haddock[0]
        pre_residue = int(args.haddock[1])
    else:
        haddock_out = args.haddock[1]
        pre_residue = int(args.haddock[0])

    #Deal with segid
    if args.segid:
        segid_1, segid_2 = args.segid
    else:
        segid_1, segid_2 = 'A','B'

    '''
    #Deal with cutoff
    if args.cutoff:
        cutoff = args.cutoff[0]
    else:
        cutoff = 0.33

    #Lastly, deal with a potential offset
    if args.offset:
        offset = args.offset[0]
    else:
        offset = 0
    '''

    with open(haddock_out, 'w') as fo:
        for peak in data_table:

            if peak[4] < args.cutoff:
                resi, resn, name = int(peak[1]), peak[2], peak[3]
                resi += args.offset

                ratio = peak[4]

                width = exp_list[np.where(peak[0] == exp_list['assign'])]['width']
                r2 = R2_from_linewidth(width)

                pre = calc_pre_from_inten_ratio(r2, peak[4])
                pre = pre / args.frac_label
                dist = calc_dist(pre, 80e-9, 600e6) #tauc and B0 hard coded to be 80 ns and 600 MHz
                dist += 7.6 #Add extra for the size of the spin label since
                            #distance restraint below is between P and methyl C and we are not
                            #modelling in the spin label

                fo.write('assign (resid {:d} and name {:s} and segid {:s}) (resid {:d} and name P and segid {:s}) {:.1f} {:.1f} {:.1f} !{:.3f}\n'.
                         format(resi, name, segid_1, pre_residue, segid_2, dist, 0.66*dist, 0.33*dist, ratio))


#If not pymol or haddock
#Make histogram of intensity ratio and print values
if not args.haddock and not args.pymol_bfactor:
    import matplotlib.pyplot as plt

    for residue in data_table:
        resi, resn, name = residue[1], residue[2], residue[3]
        int_ratio = residue[4]
        print('{:4s} {:3s} {:3s} {:.3f}'.format(resi, resn, name, int_ratio))

    dt = np.dtype([('a','<U32'), ('b','i1'), ('c', '<U32'), ('d','<U32'), ('e','f4')])
    data_table = np.asarray(data_table, dtype=object)

    if args.cutoff:
        cutoff = args.cutoff[0]
        cutoff_table = data_table[data_table[:,-1] < cutoff]

        print(cutoff_table[np.argsort(-cutoff_table[:,4])])

    else:
        print(data_table[np.argsort(-data_table[:,4])])

    n, bins, patches = plt.hist(data_table[:,4], 20, alpha = 0.75)
    #n, bins, patches = plt.hist(data_table[:,4].astype(np.float), 20, alpha = 0.75)
    plt.xlabel('Intensity Ratio')
    plt.gca().invert_xaxis()
    plt.ylabel('Counts')
    plt.grid(True)
    plt.show()
