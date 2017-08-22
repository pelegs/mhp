# Calculates mhp values for atoms in a protein,
# puts the result in the beta column of the pdb file
# and saves it as a new file.
# Written by Peleg Bar Sapir for AG Morginsky, TU-Berlin
# To do: add ability to read/write frames, neighbor cells

import sys
from pylab import *
from prody import *
import argparse
import itertools
import mhplib

def point_sphere (c=np.array([.0,.0,.0]), r=1.0, N=5):
    return [np.array([r*np.sin(t)*np.cos(phi),r*np.sin(t)*np.sin(phi),r*np.cos(t)] + c)
            for t in np.arange(.0,2*np.pi,2*np.pi/N)
            for phi in np.arange(.0,np.pi,np.pi/N)]

def mhp (p1=np.array([.0,.0,.0]), p2=np.array([.0,.0,.0]), alpha=.5, f=.0, max_r=5.0):
    max_r = 10.0 * alpha
    r = np.linalg.norm(p2-p1)
    if r <= max_r:
        return f * np.exp(-1.0 * alpha * r)
    else: return .0

parser = argparse.ArgumentParser (description="Calculates MHP for protein surface")
parser.add_argument ('-i','--input', help='Input protein name (without extension)', required=True)
parser.add_argument ('-d','--subdir', help='Input sub-directory name', required=False, default='')
parser.add_argument ('-s','--select', help='Select only part of the molecule (e.g. protein, residue, water)', required=False, default='protein')
parser.add_argument ('-p','--points', help='Number of points per atom for calculation', required=False, default=10)
parser.parse_args()
args = vars (parser.parse_args())

input_file = args['input']
subdir = args['subdir']
select = args['select']
N_points = args['points']

pdb = parsePDB(subdir + '/' + input_file + '.pdb')
psf = parsePSF(subdir + '/' + input_file + '.psf')
pdb_selection = pdb.select(select)
psf_selection = psf.select(select)
print 'Number of atoms in sub-selection: ', len(pdb_selection)

molecule = [[a.getCoords(), mhplib.vdw_radii[a.getElement()], mhplib.F_val[b.getType()]] for a, b in zip(pdb_selection, psf_selection)]
A = molecule[0]
mhp_values = [sum ([mhp(A[0], B[0], .5, B[2]) for B in molecule if B is not A]) for A in molecule]
for atom, mhp_val in zip(pdb_selection, mhp_values):
    atom.setBeta (mhp_val)

writePDB (subdir + '/' + input_file + '_mhp.pdb', pdb_selection)
