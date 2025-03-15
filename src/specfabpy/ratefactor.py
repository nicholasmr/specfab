#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2024

"""
Crystal process rate factors
"""

import numpy as np
import copy, sys, os, code # code.interact(local=locals())
from .specfabpy import specfabpy as sf__ # sf private copy 
lm__, nlm_len__ = sf__.init(20)

"""
Ice
"""

def Gamma0_exp(ugrad, T, A, Q):
    """
    DDRX rate factor as temperature-activation function (Lilien et al., 2023)
    """
    R = 8.314 # gas constant (J/mol*K)
    D = (ugrad + ugrad.T)/2
    epsE = np.sqrt(np.einsum('ij,ji',D,D)/2)
    Gamma0 = epsE * A * np.exp(-Q/(R*T))
    return Gamma0
    
def Gamma0_Lilien23_lab(ugrad, T, T0=273.15):
    """
    DDRX rate factor from lab experiments (Lilien et al., 2023)
    """
    A = 1.91e7
    Q = 3.36e4
    return Gamma0_exp(ugrad, T+T0, A, Q)

def Gamma0_Lilien23_EDC(ugrad, T, T0=273.15):
    """
    DDRX rate factor from EDC ice-core calibration (Lilien et al., 2023)
    """
    A = 4.3e7
    Q = 3.36e4
    return Gamma0_exp(ugrad, T+T0, A, Q)

"""
Olivine
"""

### TBD

