import numpy as np
cimport numpy as np
from libc.math cimport exp, sqrt 
#cython: boundscheck=False, wraparound=False, nonecheck=False

# NOTE: ndim refers to array dim, not vector dim!
# Hence an NxN array has ndim=2
# and an NxNxN array has ndim=3

# Creates N points in a sphere around
# the center, with radius r
def points_sphere(np.ndarray[float, ndim=1] centre,
                  double radius,
                  int N):
    cdef np.ndarray[double, ndim=2] sphere = np.array([
                                        [np.sin(t)*np.cos(f),
                                         np.sin(t)*np.sin(f),
                                         np.cos(t)]
                                        for t in np.arange(.0, 2*np.pi, 2*np.pi/N)
                                        for f in np.arange(.0, np.pi, np.pi/N)
                                        ])
    return centre + radius * sphere

# Actual MHP calculation between two points p1, p2
def mhp (np.ndarray[double, ndim=1] p1,
         np.ndarray[float, ndim=1] p2, 
         double f_i,
         double alpha):
    
    cdef double r = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)
    
    return f_i * exp(-alpha * r)
