""" Calculates mhp values for atoms in a protein,
    puts the result in the beta column of the pdb file
    and saves it as a new file.
    Written by Peleg Bar Sapir for AG Morginsky, TU-Berlin
    TODO: probe size seems problematic
"""

import sys
import os
from prody import *
import argparse
import itertools
import mhplib
import mhplib_c
import numpy as np

parser = argparse.ArgumentParser (description="Calculates MHP for protein surface")
parser.add_argument('-dcd','--dcd_file', help='Input dcd file', required=False)
parser.add_argument('-pdb','--pdb_file', help='Input pdb file', required=False)
parser.add_argument('-psf','--psf_file', help='Input psf file', required=True)
parser.add_argument('-select','--subselect', help='Select only part of the molecule (e.g. protein, residue, water)', required=False, default='protein')
parser.add_argument('-points','--points', help='Number of points per atom for calculation', required=False, default=64)
parser.add_argument('-probe','--probe', help='Size of probe for SAS calculations', required=False, default=1.4)
parser.add_argument('-cutoff','--cutoff', help='Cutoff distance', required=False, default=4)
parser.add_argument('-frames','--frames', help='Range of frames to use', required=False, default='1,1')
parser.parse_args()
args = vars (parser.parse_args())

num_points = int(args['points'])
cutoff_dist = float(args['cutoff'])
probe = float(args['probe'])
psf_file = args['psf_file']
pdb_file = args['pdb_file']
dcd_file = args['dcd_file']
files = [a+' ' if a else '' for a in [psf_file, pdb_file, dcd_file]]
frame_list = [int(f) for f in args['frames'].split(',')]
num_frames = frame_list[-1] - frame_list[0]
if len(frame_list) is not 2:
    print('Frame range should be 2: first frame and last frame.')
    exit(1)

print('Number of points = ', num_points)
print('Cutoff distance = ', cutoff_dist)
print('Selection is \'', args['subselect'], '\'')
print('Calculating for frames ', frame_list[0], 'to', frame_list[-1])
print('Parsing files:', ''.join(files))

# PSF variables
selection = args['subselect']
psf = parsePSF(psf_file)
selected_psf = psf.select(selection)

# Deciding between single frame (pdb) and multiple frames (dcd)
if args['pdb_file']:
    pdb = parsePDB(pdb_file)
    selected_pdb = pdb.select(selection)

    output_file = mhplib.format(pdb_file) \
                  + ' ' \
                  + mhplib.format(selection) \
                  + ' points=' \
                  + str(num_points) \
                  + ' cutoff=' \
                  + str(cutoff_dist) \
                  + ' probe=' \
                  + str(probe) \
                  + '.pdb'
    print('Number of atoms in sub-selection', selection, 'is', len(selected_psf))

    coords = selected_pdb.getCoords()
    f_vals = [mhplib.F_val[typ] for typ in selected_psf.getTypes()]
    radii  = [mhplib.vdw_radii[element[0]]+probe for element in selected_psf.getTypes()]
    molecule = [{'coords':c, 'f_val':f, 'radius':r}
                for c, f, r in zip(coords, f_vals, radii)]
    mhp_vals, total_SAS_area, positive_MHP, negative_MHP = mhplib_c.MHP_mol(molecule, coords, cutoff_dist, num_points, probe)
    print('Estimated values (log P, log P normalized, sum positive MHP values, sum negative MHP values:',
          np.sum(mhp_vals),
          np.sum(mhp_vals)/total_SAS_area,
          positive_MHP,
          negative_MHP)
    writePDB(output_file, atoms=selected_pdb, beta=mhp_vals)

if args['dcd_file']:
    traj = Trajectory(dcd_file)
    traj.link(psf)
    selected_atoms = psf.select(selection)
    traj.setAtoms(selected_atoms)

    output_files = ' '.join(['temp{}.pdb'.format(i)
                             for i in range(frame_list[0], frame_list[1]+1)])
    print('Number of atoms in sub-selection', selection, 'is', len(selected_atoms))

    logP, normalized_logP, logP_all_points, positive_vals, negative_vals = [], [], [], [], []
    for i, frame in enumerate(traj, start=frame_list[0]):
        if i > frame_list[-1]:
            break
        print('Frame', i, 'of', frame_list[-1])
        coords = np.array(frame.getAtoms().getCoords(), dtype=np.float64)
        f_vals = np.array([mhplib.F_val[typ] for typ in frame.getAtoms().getTypes()], dtype=np.float64)
        radii = np.array([mhplib.vdw_radii[element[0]]+probe for element in frame.getAtoms().getTypes()], dtype=np.float64)
        molecule = [{'coords':c, 'f_val':f, 'radius':r}
                    for c, f, r in zip(coords, f_vals, radii)]

        mhp_vals, total_SAS_area, all_mhp_vals, positive_MHP, negative_MHP = mhplib_c.MHP_mol(molecule, coords, cutoff_dist, num_points, probe)

        writePDB('temp{}.pdb'.format(i), atoms=frame.getAtoms(), beta=mhp_vals)

        logP.append(sum(mhp_vals))
        normalized_logP.append(sum(mhp_vals)/total_SAS_area)
        logP_all_points.append(sum(all_mhp_vals)/total_SAS_area)
        positive_vals.append(positive_MHP)
        negative_vals.append(negative_MHP)

    print('Estimated values (log P, log P normalized, sum positive MHP values, sum negative MHP values:\n',
          np.average(mhp_vals),        np.std(mhp_vals),
          np.average(normalized_logP), np.std(normalized_logP),
          np.average(logP_all_points), np.std(logP_all_points),
          np.average(positive_vals),   np.std(positive_vals),
          np.average(negative_vals),   np.std(negative_vals))

    print('Creating one pdb file...')
    awk_cmd = "awk 'FNR==1 && NR!=1 {print \"END\"}{print}' " \
            + output_files \
            + ' > ' \
            + mhplib.format(dcd_file) \
            + '\ ' \
            + mhplib.format(selection) \
            + '\ frames=' \
            + str(frame_list[0]) \
            + '_' \
            + str(frame_list[-1]) \
            + '\ points=' \
            + str(num_points) \
            + '\ cutoff=' \
            + str(cutoff_dist) \
            + '\ probe=' \
            + str(probe) \
            + '.pdb'
    print(awk_cmd)
    os.system(awk_cmd)
    print('Deleteing temp files...')
    os.system('rm ' + output_files)

print('')
print('Finished MHP calculation.')
