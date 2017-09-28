# Calculates mhp values for atoms in a protein,
# puts the result in the beta column of the pdb file
# and saves it as a new file.
# Written by Peleg Bar Sapir for AG Morginsky, TU-Berlin
# To do: change writing back method to not be horrible ;)

import os
import sys
from prody import *
import argparse
import itertools
import mhplib
import numpy as np
from itertools import chain

def point_sphere (c=np.array([.0,.0,.0]), r=1.0, N=5):
    return [np.array([r*np.sin(t)*np.cos(phi),r*np.sin(t)*np.sin(phi),r*np.cos(t)] + c)
            for t in np.arange(.0,2*np.pi,2*np.pi/N)
            for phi in np.arange(.0,np.pi,np.pi/N)]

def neighbor_cells(cells, Ns):
    x, y, z = cells[0], cells[1], cells[2]
    Nx, Ny, Nz = Ns[0], Ns[1], Ns[2]
    return [(i, j, k)
            for i in range(x-1, x+2) if 0 <= i < Nx
            for j in range(y-1, y+2) if 0 <= j < Ny 
            for k in range(z-1, z+2) if 0 <= k < Nz ]

def mhp (p1=np.zeros(3), p2=np.zeros(3), f=.0, alpha=.5):
    r = np.linalg.norm(p2-p1)
    return f * np.exp(-1.0 * alpha * r)

parser = argparse.ArgumentParser (description="Calculates MHP for protein surface")
parser.add_argument ('-pdb','--pdb_file', help='Input pdb file', required=True)
parser.add_argument ('-psf','--psf_file', help='Input psf file', required=True)
parser.add_argument ('-d','--dir',      help='Input directory of files', required=False, default='')
parser.add_argument ('-s','--select', help='Select only part of the molecule (e.g. protein, residue, water)', required=False, default='protein')
parser.add_argument ('-p','--points', help='Number of points per atom for calculation', required=False, default=10)
parser.add_argument ('-f','--frames', help='Frame range for calculation', required=False, type=str)
parser.parse_args()
args = vars (parser.parse_args())

pdb_file = args['pdb_file']
psf_file = args['psf_file']
directory = args['dir']
selection = args['select']
N_points = int(args['points'])
inv_N = 1.0/N_points
num_cells = 3*[50]

pdb = parsePDB(directory + '/' + pdb_file)
psf = parsePSF(directory + '/' + psf_file)
selected_pdb = pdb.select(selection)
selected_psf = psf.select(selection)
print('Files loaded and pasred.')
print('Number of atoms in sub-selection', selection, 'is', len(selected_pdb), len(selected_psf))

coords = selected_pdb.getCoords()
f_vals = [mhplib.F_val[typ] for typ in selected_psf.getTypes()]
radii  = [mhplib.vdw_radii[element[0]] for element in selected_psf.getTypes()]
molecule = [{'coords':x, 'f_val':y, 'radius':z} for x, y, z in zip(coords, f_vals, radii)]

mins = np.array([np.min(coords[:,[i]]) for i in range(3)])
maxs = np.array([np.max(coords[:,[i]]) for i in range(3)])
lengths = np.array([x1-x0 for x0, x1 in zip(mins, maxs)])
cells = [[[[] 
         for _ in range(num_cells[0]+1)]
         for _ in range(num_cells[1]+1)]
         for _ in range(num_cells[2]+1)]
for atom in molecule:
    atom['cell'] = [ int(np.floor(N*(x-m)/L)) for x, m, L, N, in zip(atom['coords'], mins, lengths, num_cells) ]
    cells[atom['cell'][0]][atom['cell'][1]][atom['cell'][2]].append(atom)
for atom in molecule:
    c = list(chain(*[ cells[i][j][k] for (i,j,k) in neighbor_cells(atom['cell'], num_cells) ]))
    atom['neighbors'] = [a for a in c if a is not atom]

mhp_list = []
for j, atom in enumerate(molecule):
    points = point_sphere(atom['coords'], atom['radius'], N_points)
    mhp_list.append( sum([ mhp(p, B['coords'], B['f_val']) for p in points for B in atom['neighbors'] ]) * inv_N )
    sys.stderr.write('\rcalculating for atom {} of {} ({} points per atom): {}   '.format(j+1, selected_pdb.numAtoms(), N_points, mhp_list[j]))

out_file = directory + '/' + pdb_file + '_' + selection + '_' + str(N_points) 
writePDB(out_file, atoms=selected_pdb, beta=mhp_list)
