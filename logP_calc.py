import numpy as np
from prody import *
import argparse
import mhplib

parser = argparse.ArgumentParser (description="Calculate log P from surface MHP (as beta values)")
parser.add_argument('-pdb','--pdb_file', help='Input PDB file', required=True)
parser.parse_args()
args = vars (parser.parse_args())

pdb_file = mhplib.rmbs(args['pdb_file'])
pdb = parsePDB(pdb_file)
print(pdb_file, sum(pdb.getBetas()))
