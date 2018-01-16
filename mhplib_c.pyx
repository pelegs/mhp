from progressbar import *
import numpy as np
cimport numpy as np
from libc.math cimport exp, sqrt 
import itertools
#cython: boundscheck=False, wraparound=False, nonecheck=False

# NOTE: ndim refers to array dim, not vector dim!
# Hence an NxN array has ndim=2
# and an NxNxN array has ndim=3

f_type = np.float64

"""
Euclidian distance
"""
cdef double sqr_distance(np.ndarray[double, ndim=1] p1,
                        np.ndarray[double, ndim=1] p2):
    return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2

"""
Using Vogel's method for generating N evenly-spaced
points on the surface of a sphere
"""
cdef np.ndarray[double, ndim=2] points_sphere(np.ndarray[double, ndim=1] center,
                                              double radius,
                                              int N):
    
    cdef double golden_angle = np.pi * (3 - np.sqrt(5))
    
    cdef np.ndarray[double, ndim=1] z = np.linspace(1 - 1.0/N, 1.0/N - 1, N).astype(f_type)
    cdef np.ndarray[double, ndim=1] r = np.sqrt(1 - z**2)
    
    cdef np.ndarray[double, ndim=1] theta = golden_angle * np.arange(N, dtype=f_type)
    cdef np.ndarray[double, ndim=1] c = np.cos(theta)
    cdef np.ndarray[double, ndim=1] s = np.sin(theta)
    
    cdef np.ndarray[double, ndim=2] points = np.zeros((N, 3), dtype=f_type)
    points[:,0] = radius * r * c
    points[:,1] = radius * r * s
    points[:,2] = radius * z
    points[:,] += center

    return points

"""
Returns only the points which lie on the atom SAS,
by the Shrake and Rupley (1973) method
"""
cdef np.ndarray[double, ndim=2] SAS(np.ndarray[double, ndim=2] points,
                                    neighbors):
    cdef np.ndarray[double, ndim=2] SAS_points = np.zeros((2,3), dtype=f_type)
    for point in points:
        keep = True
        for atom in neighbors:
            if sqr_distance(point.astype(f_type), atom['coords'].astype(f_type)) <= atom['radius']**2:
                keep = False
                break
        if keep:
            SAS_points = np.vstack((SAS_points, point)).astype(f_type)
    SAS_points = np.delete(SAS_points, (0,1), axis=0)
    return SAS_points

"""
Actual MHP calculation between two points p1, p2
"""
cdef double mhp (np.ndarray[double, ndim=1] p1,
                 np.ndarray[double, ndim=1] p2, 
                 double f_i,
                 double alpha):
    
    cdef double r = sqrt(sqr_distance(p1, p2)) 
    return f_i * exp(-alpha * r)

def neighbor_cells(indx, Ns):
    x, y, z = indx[0], indx[1], indx[2]
    Nx, Ny, Nz = Ns[0], Ns[1], Ns[2]
    return [[i, j, k]
            for i in range(x-1, x+2) if 0 <= i < Nx
            for j in range(y-1, y+2) if 0 <= j < Ny
            for k in range(z-1, z+2) if 0 <= k < Nz ]

def MHP_mol(molecule, coords, cutoff_dist, num_points, probe):
    num_atoms = len(molecule)
    mins = np.array([np.min(coords[:,[i]]) for i in range(3)])
    maxs = np.array([np.max(coords[:,[i]]) for i in range(3)])
    lengths = np.array([x1-x0 for x0, x1 in zip(mins, maxs)])
    num_cells = [int(np.ceil(L/cutoff_dist)) for L in lengths]
    
    cells = [[[[]
             for _ in range(num_cells[2]+1)]
             for _ in range(num_cells[1]+1)]
             for _ in range(num_cells[0]+1)]
    
    for atom in molecule:
        atom['cell'] = [ int(np.floor(N*(x-m)/L))
                         for x, m, L, N, in
                         zip(atom['coords'], mins, lengths, num_cells) ]
        a, b, c = atom['cell']
        cells[a][b][c].append(atom)
    
    for atom in molecule:
        neighbors = [ cells[i][j][k] for (i,j,k) in neighbor_cells(atom['cell'], num_cells) ]
        neighbor_list = list(itertools.chain(*neighbors))
        atom['neighbors'] = [a for a in neighbor_list if a is not atom]

    cdef np.ndarray[double, ndim=1] mhp_vals = np.zeros(num_atoms, dtype=f_type)
    cdef np.ndarray[double, ndim=2] points
    cdef np.ndarray[double, ndim=2] SAS_points
    cdef double atomic_SAS_area
    cdef double total_SAS_area = 0.0
    cdef double positive_MHP = 0.0
    cdef double negative_MHP = 0.0
    all_mhp_vals = []

    bar = ProgressBar(max_value=num_atoms)
    for j, atom in enumerate(molecule):
        points = points_sphere(atom['coords'].astype(f_type),
                               atom['radius'] + probe,
                               num_points)
        SAS_points = SAS(points, atom['neighbors'])
        atomic_SAS_area = 4 * np.pi * atom['radius']**2 * len(SAS_points) / num_points
        if len(SAS_points) > 0:
            area_per_point = atomic_SAS_area / len(SAS_points)
        else:
            area_per_point = 0.0
        total_SAS_area += atomic_SAS_area
        atomic_mhp = [ mhp(p, B['coords'], B['f_val'], 0.5)
                       for p in SAS_points
                       for B in atom['neighbors'] ]
        for p in atomic_mhp:
            all_mhp_vals.append(p * area_per_point)

        mhp_vals[j] = atomic_SAS_area * np.average(atomic_mhp)
        if np.isnan(mhp_vals[j]):
            mhp_vals[j] = 0.0
        if mhp_vals[j] > 0.0:
            positive_MHP += mhp_vals[j]
        if mhp_vals[j] < 0.0:
            negative_MHP += mhp_vals[j]

        bar.update(j)
    
    print('')

    return mhp_vals, total_SAS_area, all_mhp_vals, positive_MHP, negative_MHP

def print_SAS_points(molecule, coords, cutoff_dist, num_points, probe):
    cdef np.ndarray[double, ndim=2] SAS_points
    cdef np.ndarray[double, ndim=2] points
    
    num_atoms = len(molecule)
    mins = np.array([np.min(coords[:,[i]]) for i in range(3)])
    maxs = np.array([np.max(coords[:,[i]]) for i in range(3)])
    lengths = np.array([x1-x0 for x0, x1 in zip(mins, maxs)])
    num_cells = [int(np.ceil(L/cutoff_dist)) for L in lengths]

    cells = [[[[]
             for _ in range(num_cells[2]+1)]
             for _ in range(num_cells[1]+1)]
             for _ in range(num_cells[0]+1)]
    
    for atom in molecule:
        atom['cell'] = [ int(np.floor(N*(x-m)/L))
                         for x, m, L, N, in
                         zip(atom['coords'], mins, lengths, num_cells) ]
        a, b, c = atom['cell']
        cells[a][b][c].append(atom)
    
    for atom in molecule:
        neighbors = [ cells[i][j][k] for (i,j,k) in neighbor_cells(atom['cell'], num_cells) ]
        neighbor_list = list(itertools.chain(*neighbors))
        atom['neighbors'] = [a for a in neighbor_list if a is not atom]

    for j, atom in enumerate(molecule, probe):
        points = points_sphere(atom['coords'].astype(f_type),
                               atom['radius'] + probe,
                               num_points)
        for point in SAS(points, atom['neighbors']):
            print('C', ' '.join(map(str, point)))
