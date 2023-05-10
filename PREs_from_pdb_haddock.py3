#!/usr/bin/env python3
#
# Calculate distances from a pdb file
#
# Uses HADDOCK format unambiguous distance file to know which atoms to
# calculate the distances
#
# MP Latham Aug. 9, 2021

import sys, re
import argparse

import numpy as np
import matplotlib.pyplot as plt

from math import sqrt
from scipy.constants import mu_0, pi

try:
    from termcolor import colored
except ImportError:
    print("\nFor color output in terminal window run")
    print("python3 -m pip install --upgrade termcolor\n")

parser = argparse.ArgumentParser('Script to calculate distances from PDB file')
parser.add_argument('-pdb', type=str, required=True, nargs='+', help='pdb file')
parser.add_argument('-tbl', type=str, required=True, help='HADDOCK format table of atom pairs and distances')

if len(sys.argv) == 1:
    parser.print_help()
    exit()

args = parser.parse_args()

#Functions
def read_pdb(pdb):
    '''Takes pdb file and returns numpy array of 
       resi, resn, atom name, x-coor, y-coor, and z-coor'''

    atoms = list()

    with open(pdb, 'r') as pdb_file:
        lines = pdb_file.readlines()

        for line in lines:
            entry = line.split()

            if entry[0] == 'TER':
                continue
            if entry[0] != 'ATOM':
                continue

            #resi, resn, atom name, chain, x, y, z
            if len(entry) == 12:
                resi = int(entry[4][1:5])
                chain = str(entry[4][0])

                atoms.append([resi, entry[3], entry[2], chain,
                    float(entry[5]), float(entry[6]), float(entry[7])])
            
            else: 
                if isinstance(entry[4], str):
                    atoms.append([int(entry[5]), entry[3], entry[2], entry[4],
                                  float(entry[6]), float(entry[7]), float(entry[8])])

    return np.array(atoms)

def calc_dist(x1, y1, z1, x2, y2, z2):
    '''Calculate the distance between two atoms'''

    dist = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    return float('{0:.2f}'.format(dist))

def calc_pre(dist, tauc, B0):
    '''Calculate the PRE between two atoms based on the distance.'''

    gamma_H = 26.7522128e+07 #1H gyromagnetic ratio (T^-1 s^-1)
    g = -2.0023 #electron g value
    beta = 9.27401e-24 #Bohr magneton (J T^-1)
    spin = 1./2. #nitroxide spin label

    C = (mu_0**2) * (gamma_H**2) * (g**2) * (beta**2) * spin*(spin + 1)
    C = C/(15 * (4 * pi)**2) #Constants
    
    w_h = 2 * pi *B0
    B = (4.0*tauc + 3.0*tauc / (1 + (tauc * w_h)**2)) #Spectral density part

    dist = dist*1e-10 #convert Angstoms to meter

    return C*(dist**(-6.))*B #PRE


#Deal with the possibility of having wild card in list of pdbs
if '*' in args.pdb:
    import glob
    pdbs = glob.glob(args.pdb)
else:
    pdbs = args.pdb

master_list = list()
#Loop over structures
for structure in pdbs:

    exp_dist = list()
    pairs = list()
    temp_list = list()

    #Read in the PDB
    pdb = read_pdb(structure)

    #Open the file containing pairs of residues
    with open(args.tbl, 'r') as distance_file:
        lines = distance_file.readlines()
    
        for line in lines:

            first, second = re.findall('\(([^)]+)', line) #get text between '()'

            resi1 = re.search('resid(\s+\w+)', first).groups()[0].strip() #get resi, which is the number after resid in file
            name1 = re.search('name(\s+\w+)', first).groups()[0].strip()  #get name, which is the string after name in file
            chain1 = re.search('segid(\s+\w+)', first).groups()[0].strip() #get chain, which is the string after segid in file 

            if pdb[np.where((pdb[:,0] == resi1) & (pdb[:,2] == name1) & (pdb[:,3] == chain1))].size > 0: #find this atom in pdb file
                ca1 = pdb[np.where((pdb[:,0] == resi1) & (pdb[:,2] == name1) & (pdb[:,3] == chain1))][0] 
                x1, y1, z1 = float(ca1[4]), float(ca1[5]), float(ca1[6]) #get x,y,z coordinates

                resi2 = re.search('resid(\s+\w+)', second).groups()[0].strip() #do the same for the second atom
                name2 = re.search('name(\s+\w+)', second).groups()[0].strip()
                chain2 = re.search('segid(\s+\w+)', second).groups()[0].strip()

                ca2 = pdb[np.where((pdb[:,0] == resi2) & (pdb[:,2] == name2) & (pdb[:,3] == chain2))][0]
                x2, y2, z2 = float(ca2[4]), float(ca2[5]), float(ca2[6])

                dist = calc_dist(x1, y1, z1, x2, y2, z2)
                '''Right now distance calculation is only happening for methyl
                carbon and not the average of the 3 methyl protons. Our structures
                do not have 1H and HADDOCK does not add them in.'''

                pre = calc_pre(dist, 80e-9, 600e6) #tauc and 1H frequency hard coded for MR complex and 600 MHz NMR 

                line_temp = line.split()

                if '!' in line_temp[-1]: #Newer version of HADDOCK generating script includes a comment of the intensity ratio or PRE value
                    e_dist = float(line_temp[-4])
                    l_bound = float(line_temp[-3])
                    u_bound = float(line_temp[-2])

                else:
                    e_dist = float(line_temp[-3])
                    l_bound = float(line_temp[-2])
                    u_bound = float(line_temp[-1])
            
                pairs.append([resi1, resi2])
                exp_dist.append([e_dist])

                max_dist = e_dist + u_bound
                min_dist = e_dist - l_bound

                #Check if termcolor has been loaded
                if 'termcolor' in sys.modules:
                    #If the distance from the PDB is within the bounds of the HADDOCK input
                    #this will print the distance out to the terminal in green.
                    if dist > min_dist and dist < max_dist:
                        colored_dist = colored(str(dist), "green")
                        temp_list.append(colored_dist)
                    else:
                        temp_list.append(dist) #Otherwise it will be black

                else:
                    temp_list.append(dist)
            else:
                pass #If residue in HADDOCK table is not in the PDB, pass.

    master_list.append(temp_list)

master_list = np.array(master_list)
#master_list = np.array(master_list, dtype=float)
pairs = np.array(pairs, dtype=int)
exp_dist = np.array(exp_dist,  dtype=float)

if len(exp_dist) != 0:
    header = 'resi_1 resi_2 exp_dist ' + ''.join([str(i)+' ' for i in pdbs])
    print('\n')
    print(header)

    for rows in np.hstack((pairs, exp_dist, master_list.T)):
        row = [str(i) for i in rows]
        print("{0}".format(' '.join(map(str, row))))
else:
    header = 'resi_1 resi_2 ' + ''.join([str(i)+' ' for i in pdbs])
    print('\n')
    print(header)

    for rows in np.hstack((pairs, master_list.T)):
        row = [str(i) for i in rows]
        print("{0}".format(' '.join(map(str, row))))


