#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>

"""
Backconstruct a4 from its eigenframe (eigentensors and eigenvalues)
"""

import copy, sys, os, code # code.interact(local=locals())
import numpy as np
from specfabpy import specfab as sf
from specfabpy import integrator as sfint

L = 8
lm, nlm_len = sf.init(L) 

def print_summary(nlm):
    nlm = nlm[-1,:] 
    a4 = sf.a4(nlm)
    *Qi, eigvals = sf.a4_eigentensors(nlm)
    a4_test = sum([eigvals[ii]*np.einsum('ij,kl->ijkl', Qi[ii], Qi[ii]) for ii in range(6)])
    
    a4mat = sf.a4_to_mat(a4)
    a4mat_test = sf.a4_to_mat(a4_test)
    
    diff = np.abs(a4mat - a4mat_test)
    print(a4mat)
    print(a4mat_test)
    print('min diff %.2e :: max diff %.2e :: sum diff %.2e'%(np.amin(diff), np.amax(diff), np.sum(diff)))

# some CPO state...
Nt = 100
kwargs_LROT = dict(iota=1, Gamma0=None, nu=1)

mod = dict(type='ss', plane=2)
strain_target = np.deg2rad(50) # simulate  parcel deformation until this target strain
nlm, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, **kwargs_LROT)
print_summary(nlm)

mod = dict(type='ps', r=-1, T=-1, axis=0)
strain_target = +4 # simulate  parcel deformation until this target strain
nlm, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, **kwargs_LROT)
print_summary(nlm)

