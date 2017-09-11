# Calculates mhp values for atoms in a protein,
# puts the result in the beta column of the pdb file
# and saves it as a new file.
# Written by Peleg Bar Sapir for AG Morginsky, TU-Berlin
# To do: change writing back method to not be horrible :P

import os
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
selection = args['select']
N_points = int(args['points'])
inv_N = 1.0/N_points
cutoff_dist = 5.0

traj = Trajectory(subdir + '/' + input_file + '.dcd')
psf = parsePSF(subdir + '/' + input_file + '.psf')
traj.link(psf)
selected_atoms = psf.select(selection)
traj.setAtoms(selected_atoms)
print 'Files loaded and pasred.'
print 'Number of atoms in sub-selection', selection, 'is', len(selected_atoms)
num_frames = len(traj)
frame_list = ' '.join(['temp{}.pdb'.format(i) for i in range(num_frames)])

for i, frame in enumerate(traj):
    print 'Frame', i, '...'
    coords = frame.getAtoms().getCoords()
    f_vals = [mhplib.F_val[typ] for typ in frame.getAtoms().getTypes()]
    radii = [mhplib.vdw_radii[element[0]] for element in frame.getAtoms().getTypes()]
    distances = buildDistMatrix(frame)
    molecule = [{'coords':x, 'f_val':y, 'radius':z} for x, y, z in zip(coords, f_vals, radii)]
    for j, atom in enumerate(molecule):
        atom['neighbors'] = [ molecule[i] for j,d in enumerate(distances[i]) if d <= cutoff_dist ]

    mhp_list = []
    for j, atom in enumerate(molecule):
        points = point_sphere(atom['coords'], atom['radius'], N_points)
        mhp_list.append( sum([ mhp(atom['coords'], B['coords'], B['f_val']) for B in atom['neighbors'] ]) * inv_N )
        sys.stderr.write('\rcalculating for atom {} of {} ({} points per atom)   '.format(j+1, len(selected_atoms), N_points))

    writePDB('temp{}.pdb'.format(i), atoms=frame.getAtoms(), beta=mhp_list)
    print ''

print 'Done.'
print 'Saving to file {}'.format(input_file)

print 'Creating one pdb file...'
os.system('awk \'FNR==1 && NR!=1 {print "END"}{print}\' ' + frame_list + ' > ' + subdir + '/' + input_file + '_mhp.pdb')
print 'Deleteing temp files...'
os.system('rm ' + frame_list)
print 'Done.'
