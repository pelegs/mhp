"""
Calculates mhp value changes over time.
MHP values should be saved as beta values in
dcd/pdb files.
Written by Peleg Bar Sapir for AG Morginsky, TU-Berlin
"""

import sys
import os
import argparse
import mhplib
from prody import *
import subprocess
import numpy as np

parser = argparse.ArgumentParser (description="Calculates MHP for protein surface")
parser.add_argument('-pdb','--pdb_file', help='Input pdb file', required=False)
parser.add_argument('-a','--atoms', help='select atoms', required=True)
parser.add_argument('-f','--frames', help='Range of frames to use', required=False)
parser.add_argument('-o','--output', help='Output file name', required=True)
parser.parse_args()
args = vars (parser.parse_args())

pdb_file = args['pdb_file']
if args['frames']:
    frame_list = [int(f) for f in args['frames'].split(',')]
    frames = range(frame_list[0], frame_list[-1]+1)
atom_ids = [int(a) for a in args['atoms'].split(',')]

output_file = open(args['output'], 'w')

## creating separate pdb files in dir 'temp'
#subprocess.run(['mkdir', 'temp'])
#subprocess.run(['csplit', pdb_file, '/END/', '{*}', '-f', 'temp/pdb'],
#               stdout=subprocess.PIPE)

# reading odb files one after the other
files = ['temp/'+f for f in os.listdir('temp')]
files.sort()
for i, frame in enumerate(files):
    pdb = parsePDB(frame)
    mhp_avg = np.average([pdb.getBetas()[i] for i in atom_ids])
    output_file.write(str(i) + ' ' + str(mhp_avg) + '\n')

output_file.close()
