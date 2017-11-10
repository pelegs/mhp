from progressbar import *
import numpy as np
cimport numpy as np
from libc.math cimport exp, sqrt 
import itertools
#cython: boundscheck=False, wraparound=False, nonecheck=False

#NOTE: ndim refers to array dim, not vector dim!
#Hence an NxN array has ndim=2
#and an NxNxN array has ndim=3

"""
Using Vogel's method for generating N evenly-spaced
points on the surface of a sphere
"""
cdef np.ndarray[double, ndim=2] points_sphere(np.ndarray[double, ndim=1] centre,
                                              double radius,
                                              int N):
    
    cdef double golden_angle = np.pi * (3 - np.sqrt(5))
    
    cdef np.ndarray[double, ndim=1] z = np.linspace(1 - 1.0/N, 1.0/N - 1, N)
    cdef np.ndarray[double, ndim=1] r = np.sqrt(1 - z**2)
    
    cdef np.ndarray[double, ndim=1] theta = golden_angle * np.arange(N)
    cdef np.ndarray[double, ndim=1] c = np.cos(theta)
    cdef np.ndarray[double, ndim=1] s = np.sin(theta)
    
    cdef np.ndarray[double, ndim=2] points = np.zeros((N, 3))
    points[:,0] = radius * r * c
    points[:,1] = radius * r * s
    points[:,2] = radius * z
    points[:,] += centre

    return points

"""
Actual MHP calculation between two points p1, p2
"""
cdef double mhp (np.ndarray[double, ndim=1] p1,
                 np.ndarray[float, ndim=1] p2, 
                 double f_i,
                 double alpha):
    
    cdef double r = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)
    
    return f_i * exp(-alpha * r)

# Returns indices of all neighbor cells
# for a cell with index indx and Ns number of cells
# (Ns is a triple)
#cdef np.ndarray[long, ndim=2] neighbor_cells(np.ndarray[long, ndim=1] indx,
#                                             np.ndarray[long, ndim=1] Ns):
#    cdef long x = indx[0]
#    cdef long y = indx[1]
#    cdef long z = indx[2]
#    
#    np.ndarray[long, ndim=2] cells = np.array([[i, j, k]
#            for i in range(x-1, x+2) if 0 <= i < Ns[0]
#            for j in range(y-1, y+2) if 0 <= j < Ns[1]
#            for k in range(z-1, z+2) if 0 <= k < Ns[2]])
#
#    return cells

def neighbor_cells(indx, Ns):
    x, y, z = indx[0], indx[1], indx[2]
    Nx, Ny, Nz = Ns[0], Ns[1], Ns[2]
    return [[i, j, k]
            for i in range(x-1, x+2) if 0 <= i < Nx
            for j in range(y-1, y+2) if 0 <= j < Ny
            for k in range(z-1, z+2) if 0 <= k < Nz ]

def MHP_mol(molecule, coords, cutoff_dist, num_points):
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

    cdef np.ndarray[double, ndim=1] mhp_vals = np.zeros(num_atoms)
    cdef np.ndarray[double, ndim=2] points
    bar = ProgressBar(max_value=num_atoms)
    
    for j, atom in enumerate(molecule):
        points = points_sphere(atom['coords'],
                               atom['radius'],
                               num_points)
        mhp_vals[j] = sum([ mhp(p, B['coords'], B['f_val'], 0.5)
                            for p in points
                            for B in atom['neighbors'] ]) / num_points
        
        bar.update(j)
    
    print('')

    return mhp_vals
