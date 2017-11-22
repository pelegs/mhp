from prody import *
import sys
import itertools
import numpy as np
from progressbar import *
import os
import mhplib_c

# Format strings to remove extensions, escape paranthesis, etc.
def format(s):
    name = os.path.splitext(s)[0]   # remove file extension
    name = name.replace('(', '\(')
    name = name.replace(')', '\)')
    name = name.replace(' ', '\ ')
    return name

def rmbs(s):
    return s.replace('\ ', ' ')

# Returns indices of all neighbor cells
# for a cell with index indx and Ns number of cells
# (Ns is a triple)
def neighbor_cells(indx, Ns):
    x, y, z = indx[0], indx[1], indx[2]
    Nx, Ny, Nz = Ns[0], Ns[1], Ns[2]
    return [[i, j, k]
            for i in range(x-1, x+2) if 0 <= i < Nx
            for j in range(y-1, y+2) if 0 <= j < Ny
            for k in range(z-1, z+2) if 0 <= k < Nz ]

# Van der waals radii in Angstrom
vdw_radii={
        'H': 1.2,
        'Li': 1.82,
        'B': 1.8,
        'C': 1.7,
        'N': 1.55,
        'O': 1.52,
        'F': 1.47,
        'Na': 2.27,
        'Mg': 1.73,
        'Si': 2.10,
        'P': 1.80,
        'S': 1.80,
        'Cl': 1.75,
        'K': 2.75,
        'Ca': 2.40,
        'Br': 1.85,
        'I': 1.98
        }

# F values by atom type
F_val={
        'HB1': 0.6666,
        'HB2': 0.6666,
        'HA': 0.7341,
        'HA1': 0.5372,      # Amide carbon H [i.e. N--CH--C(=O)]
        'HA2': 0.5372,      # Same as previos
        'HA3': 0.5372,      # Same as previos
        'CT2A': -0.6805,    # Carbon in CHR2X
        'CT3': -0.27,
        'C': -0.2405,
        'CY': 0.2952,       # Aromatic C with no H
        'CPT': 0.2952,      # Same as previos
        'CAI': -0.3962,     # Aromatic C with 1 H
        'CD': -0.3962,      # Same as previos
        'O': -0.0233,
        'H': -0.1036,
        'CT2': -1.0120,
        'HB': 0.5234,
        'NH3': -1.4439,
        'NY': 0.1259,       # N in aromatic ring with one H attached to it
        'HP': 0.6301,
        'CT1': -0.6681,
        'HC': -0.1036,
        'NC2': -1.4439,
        'CP1': -0.6805,
        'CP2': -1.0120,
        'CP3': -1.2486,
        'CC': 0.0200,
        'NH1': 0.3168,
        'NH2': 0.5113,
        'OB': -0.0233,      # Amide O
        'OC': 0.0127,
        'N': 0.1624,
        'CA': 0.3251,
        'HS': -0.1036,
        'CTL2': -1.0120,
        'OH1': -0.0127,     # Carboxyl O
        'HAL2': 0.6666,     # H attached to SP3 carbon without X on next C
        'OHL': .0,          # Need to understant this.
        'HOL': .0,          # ...and this too.
        'CTL1': -0.6805,
        'PL': -0.1726,      # Phosphate ether P
        'O2L': -0.7941,     # negativly charged O
        'OSL': 0.0324,      # X--O--R
        'OSLP': 0.0324,     # X--O--R
        'HAL1': 0.6301,     # H in X--CR2--H
        'HEL1': 0.6301,     # H in R-CH=R
        'CEL1': -0.3962,    # C in =CHR
        'HAL3': 0.7341,     # H in R-CH3 (methyl H)
        'CTL3': -1.5603,    # C in R-CH3 (methyl C)
        'NH3L': -1.4439,    # N+
        'CL': -0.3858,      # carbon in RHCX2
        'OBL': -0.0233,     # =O
        'HCL': -0.1036,     # H in H-X
        'S': 0.6146,        # S in RSH
        'HS': -0.1036,      # H in RSH (H attached to heteroatom)
        'NR1': -0.3168,     # Aromatic N with H
        'NR2': 0.0132,      # Aromatic N with no H
        'CPH1': -0.0909,    # Aromaric C with no H but X and R
        'CPH2': -0.0244,    # Aromatic C with X and H
        'HR1': 0.5234,      # Alpha H
        'HR3': 0.5234,      # Samea as previous
        }
