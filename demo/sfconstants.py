# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

# Material constants

# from sfconstants import *

class sfconst():

    # Effective grain parameters for polycrystalline ice
    ice = {
        'viscoplastic': {
            # These are the linear mixed Taylor--Sachs best-fit parameters from Rathmann and Lilien (2021)
            'linear':    ((1,1e3), 0.0125, 1), # (Eij_grain, alpha, n_grain)
            # These are the nonlinear Sachs-only best-fit parameters (Rathmann et al., 2021) 
            'nonlinear': ((1,1e4), 0.0, 3),    # (Eij_grain, alpha, n_grain)
        },
        'elastic': {
            'linear': {},
        }    
    }

    # Effective grain parameters for polycrystalline olivine
    olivine = {
        'viscoplastic': {
            'linear':    ((1,1,1, 1,1,1e1), 0.3, 1), # (Eij_grain, alpha, n_grain)
            'nonlinear': None,
        },
        'elastic': {
            'linear': {},
        }    
    }
