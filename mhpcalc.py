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
import time

parser = argparse.ArgumentParser (description="Calculates MHP for protein surface")
parser.add_argument('-dcd','--dcd_file', help='Input dcd file', required=False)
parser.add_argument('-pdb','--pdb_file', help='Input pdb file', required=False)
parser.add_argument('-psf','--psf_file', help='Input psf file', required=True)
parser.add_argument('-s','--subselect', help='Select only part of the molecule (e.g. protein, residue, water)', required=False, default='protein')
parser.add_argument('-p','--points', help='Number of points per atom for calculation', required=False, default=10)
parser.add_argument('-c','--cutoff', help='Cutoff distance', required=False, default=4)
parser.add_argument('-f','--frames', help='Range of frames to use', required=False)
parser.parse_args()
args = vars (parser.parse_args())

num_points = int(args['points'])
cutoff_dist = float(args['cutoff'])

# PSF variables
psf_file = args['psf_file']
selection = args['subselect']
psf = parsePSF(psf_file)

# Deciding between single frame (pdb) and multiple frames (dcd)
if args['pdb_file']:
    pdb_file = args['pdb_file']
    pdb = parsePDB(pdb_file)
    selected_pdb = pdb.select(selection)
    
    output_file = mhplib.rm_ext(pdb_file) \
                  + '_' \
                  + selection \
                  + '_p' \
                  + str(num_points) \
                  + '_cutoff' \
                  + str(cutoff_dist) \
                  + '.pdb'
    print('Number of atoms in sub-selection', selection, 'is', len(selected_pdb))
    
    start_time = time.time()
    
    coords = selected_pdb.getCoords()
    f_vals = [mhplib.F_val[typ] for typ in selected_psf.getTypes()]
    radii  = [mhplib.vdw_radii[element[0]] for element in selected_psf.getTypes()]
    molecule = [{'coords':c, 'f_val':f, 'radius':r}
                for c, f, r in zip(coords, f_vals, radii)]
    mhp_vals = mhplib.MHP_mol(molecule, coords, cutoff_dist, num_points)
    writePDB(output_file, atoms=selected_pdb, beta=mhp_vals)

if args['dcd_file']:
    dcd_file = args['dcd_file']
    traj = Trajectory(dcd_file)
    traj.link(psf)
    selected_atoms = psf.select(selection)
    traj.setAtoms(selected_atoms)
    
    output_files = ' '.join(['temp{}.pdb'.format(i) for i in frame_list])
    frames = args['frames']
    frame_list = [int(frame) for frame in args['frames'].split(',')]
    if len(frame_list) is not 2:
        print('Frame range should be 2: first frame and last frame.')
        exit()
    num_frames = frame_list[-1] - frame_list[0]
    print('Number of atoms in sub-selection', selection, 'is', len(selected_atoms))

    for i, frame in enumerate(traj, start=frame_list[0]):
        if i > frame_list[-1]:
            break
        print('Frame', i, 'of', frame_list[-1])
        coords = frame.getAtoms().getCoords()
        f_vals = [mhplib.F_val[typ] for typ in frame.getAtoms().getTypes()]
        radii = [mhplib.vdw_radii[element[0]] for element in frame.getAtoms().getTypes()]
        molecule = [{'coords':c, 'f_val':f, 'radius':r}
                    for c, f, r in zip(coords, f_vals, radii)]
        mhp_vals = mhplib.MHP_mol(molecule, coords, cutoff_dist, num_points)
        writePDB('temp{}.pdb'.format(i), atoms=frame.getAtoms(), beta=mhp_vals)

    print('Creating one pdb file...')
    awk_cmd = "awk 'FNR==1 && NR!=1 {print \"END\"}{print}' " \
            + output_files \
            + ' > ' \
            + mhplib.rm_ext(dcd_file) \
            + '_' \
            + selction \
            + '_frames_' \
            + str(frame_list[0]) \
            + '_' \
            + str(frame_list[-1]) \
            + '_p' \
            + str(num_points) \
            + '_cutoff' \
            + str(cutoff_dist) \
            + '.pdb'
    print(awk_cmd)
    os.system(awk_cmd)
    print('Deleteing temp files...')
    os.system('rm ' + output_files)

elapsed_time = time.time() - start_time
print('')
print('Finished MHP calculation.')
print('Time elapsed: %02f seconds' % elapsed_time)
