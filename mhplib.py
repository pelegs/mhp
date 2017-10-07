# Creates N points in a sphere around
# the center, with radius r
def points_sphere(centre, radius, N):
    sphere = np.array([[np.sin(t)*np.cos(f),
                        np.sin(t)*np.sin(f),
                        np.cos(t)]
                        for t in np.arange(.0, 2*np.pi, 2*np.pi/N)
                        for f in np.arange(.0, np.pi, np.pi/N)])
    return centre + radius * sphere

# Actual MHP calculation between two points p1, p2
# Can choose form for distance function, default: exp(-d/2)
def mhp (p1, p2, f_i, form='Broto', alpha=.5):
    r = np.linalg.norm(p2-p1)

    if form == 'Audry':
        D = alpha / (alpha + r)
    elif form == 'Fauchere':
        D = exp(-r)
    else:
        D = exp(-alpha*r)
    
    return f_i * D

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
        'HA': 0.7341,
        'CT3': -0.27,
        'C': -0.2405,
        'O': -0.0233,
        'H': -0.1036,
        'CT2': -1.0120,
        'HB': 0.5234,
        'NH3': -1.4439,
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
        'OC': 0.0127,
        'N': 0.1624,
        'CA': 0.3251,
        'HS': -0.1036,
        'CTL2': -1.0120,
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
        'HCL': -0.1036      # H in H-X
        }
