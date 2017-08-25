# Calculates mhp values for atoms in a protein,
# puts the result in the beta column of the pdb file
# and saves it as a new file.
# Written by Peleg Bar Sapir for AG Morginsky, TU-Berlin
# To do: add ability to read/write frames

import sys
from prody import *
import argparse
import itertools
import mhplib
import numpy as np

def point_sphere (c=np.array([.0,.0,.0]), r=1.0, N=5):
    return [np.array([r*np.sin(t)*np.cos(phi),r*np.sin(t)*np.sin(phi),r*np.cos(t)] + c)
            for t in np.arange(.0,2*np.pi,2*np.pi/N)
            for phi in np.arange(.0,np.pi,np.pi/N)]

def mhp (p1=np.array([.0,.0,.0]), p2=np.array([.0,.0,.0]), f=.0, alpha=.5):
    r = np.linalg.norm(p2-p1)
    return f * np.exp(-1.0 * alpha * r)

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
N_points = int(args['points'])
inv_N = 1.0/N_points

pdb = parsePDB(subdir + '/' + input_file + '.pdb')
psf = parsePSF(subdir + '/' + input_file + '.psf')
pdb_selection = pdb.select(select)
psf_selection = psf.select(select)
print 'Number of atoms in sub-selection: ', len(pdb_selection)
print 'Building distance matrix...'
dist = buildDistMatrix(pdb_selection)
print 'Done.'

molecule = [atom for atom in pdb_selection]
cutoff_dist = 5.0
neighbors_pdb = [[pdb[pdb_selection.getIndices()[i]] for i,v in enumerate(atom) if v <= cutoff_dist] for atom in dist]
neighbors_psf = [[psf[psf_selection.getIndices()[i]] for i,v in enumerate(atom) if v <= cutoff_dist] for atom in dist]
mhp_values = []
for i, atom in enumerate(molecule):
    points = point_sphere(atom.getCoords(), mhplib.vdw_radii[atom.getElement()], N_points)
    mhp_atom = sum([mhp(p, B.getCoords(), mhplib.F_val[C.getType()]) for p in points for B,C in zip(neighbors_pdb[i], neighbors_psf[i])]) * inv_N
    mhp_values.append(mhp_atom)
    sys.stderr.write('\rcalculating for atom {} of {} ({} points per atom): {}   '.format(i, len(pdb_selection), N_points, mhp_atom))
print '\n'

for atom, mhp_val in zip(pdb_selection, mhp_values):
    atom.setBeta (mhp_val)

print 'Writing to ' + input_file + '_mhp_N_' + str(N_points) + '.pdb...'
writePDB (subdir + '/' + input_file + '_mhp_N_' + str(N_points) + '.pdb', pdb_selection)
print 'Done.'
