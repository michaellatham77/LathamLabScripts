#!/usr/bin/env python3
#
# Script to extract/reorder columns of tables a la getCols from nmrPipe
#
# MP Latham, Nov. 14, 2022
######################################

import sys
import argparse

import numpy as np

parser = argparse.ArgumentParser(description='Script extract/reorder columns from table.')
parser.add_argument('-input', required=True, type=str, nargs=1, help='Input table')
parser.add_argument('-output', type=str, nargs=1, help='Output table')
parser.add_argument('-cols', required=True, type=int, nargs="+", help='Columns to extract to reorder. The first column is "0"')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

in_file = np.genfromtxt(args.input[0], skip_header=1, usecols=args.cols, dtype=None, encoding=None)

np.savetxt(args.output[0], in_file, fmt="%s", delimiter="\t")
