#!/usr/bin/env python
# 07/12/21 - Stephan Azatian
# Changed python command to local environment
# Old:
# #!/usr/bin/python

import sys
import argparse
import numpy as np

#Argument parser
parser = argparse.ArgumentParser(prog='CSPs_deltadelta.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''\
    Script to take the difference of chemical shift perturbations
    from two tables formatted for data_2bfactor pymol script
    (i.e., CSPs made with amide_CSPs.py or methyl_CSPs.py) 

    Assumes data is in the format
    129  Val   CG1  0.0146
    ''')
parser.add_argument('-CSP1', type=str, required=True, nargs=1, help='Table of first CSPs')
parser.add_argument('-CSP2', type=str, required=True, nargs=1, help='Table of second CSPs')
parser.add_argument('-out', type=str, nargs=1, default='CSP_deltadelta.tab', help='Output file')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

#Read in the data and put into Scipy record array
csp_1 = np.genfromtxt(args.CSP1[0], comments='#', dtype=None,
            names=['resi','resn','name','csp'])
csp_2 = np.genfromtxt(args.CSP2[0], comments='#', dtype=None,
            names=['resi','resn','name','csp'])

if type(args.out).__name__ == 'list':
    outfile = open(args.out[0], 'w')
elif type(args.out).__name__ == 'str':
    outfile = open(args.out, 'w')
else:
    print('unrecognized type {:s}'.format(args.out))
    sys.exit(1)

for resi1, resn1, name1, csp1 in csp_1:
    for resi2, resn2, name2, csp2 in csp_2:

        if resi1 == resi2 and resn1 == resn2 and name1 == name2:
            print('{:4d} {:3s} {:3s} {:.3f}'.format(resi1, resn1, name1, csp1 - csp2))
            outfile.write('{:4d} {:3s} {:3s} {:.3f}\n'.format(resi1, resn1, name1, csp1 - csp2))

outfile.close()
