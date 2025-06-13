# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

"""
Effective grain parameters for polycrystalline materials
"""

ice = {
    'viscoplastic': {
        # (Eij_grain, alpha, n_grain)
        'linear':    ((1,1e3), 0.0125, 1), # Linear mixed Taylor--Sachs best-fit parameters from Rathmann and Lilien (2021)
        'linearL10': ((1,1e3), 0.455,  1), # Linear mixed Taylor--Sachs best-fit parameters from Rathmann et al. (2025)
        'nonlinear': ((1,1e4), 0.0, 3),    # Nonlinear Sachs-only best-fit parameters from Rathmann et al. (2021)
        'CAFFE':       (0.1,  10, 1), # Emin, Emax, n_grain
        'CAFFEdagger': (0.2, 3.2, 1), 
        'ESTAR':       (3, 8), # Emin, Emax
    },
    'elastic': {
        # (C11,C33,C55,C12,C13)
        'Bennett1968' : (14.060e9, 15.240e9, 3.060e9, 7.150e9, 5.880e9), \
        'Gammon1983'  : (13.929e9, 15.010e9, 3.014e9, 7.082e9, 5.765e9), \
        'Green1956'   : (13.33e9, 14.28e9, 3.26e9, 6.03e9, 5.08e9), \
        'Gagnon1988'  : (14.48e9, 15.63e9, 3.15e9, 7.35e9, 5.97e9), \
        'Dantl1968'   : (13.21e9, 14.43e9, 2.89e9, 6.69e9, 5.84e9), \
        'Bass1957'    : (13.3e9, 14.2e9, 3.06e9, 6.3e9, 4.6e9), \
        'Brockamp1964': (13.63e9, 14.85e9, 3.04e9, 6.69e9, 5.15e9), \
    },
    'electromagnetic': {
        # (epsr_para, epsr_perp, mu); relative permitivities and permeability, where para is parallel to symmetry (c) axis, and perp the orthogonal plane of isotropy.
        'Fujita2000': (1,1,1) # to be updated
    },
    'density': 917 # kg/m^3
}

olivine = {
    'viscoplastic': {
        # (Eij_grain, alpha, n_grain)
#        'linear':    ((1,1,1, 1,1,1e1), 0.3, 1), # Linear mixed Taylor--Sachs best-fit parameters
        'linear':    ((1,1,1, 1,1,1e2), 0, 1), # Linear Sachs parameters from Rathmann et al. (2023)
        'nonlinear': None,
    },
    'elastic': {
        # (C11,C22,C33,C44,C55,C66,C23,C13,C12)
        'Abramson1997': (320.5e9, 196.5e9, 233.5e9,  64.0e9, 77.0e9, 78.7e9,  76.8e9, 71.6e9, 68.15e9),
        'Jacobsen2008': (320.2e9, 195.9e9, 233.8e9,  63.5e9, 76.9e9, 78.1e9,  78.5e9, 70.5e9, 67.90e9),  
    },
    'density': 3355 # kg/m^3
}

