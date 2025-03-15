#!/usr/bin/python3
# Nicholas Rathmann

"""
Test how fabric--stress misalignment affects coaxiality
"""

import numpy as np
from numpy import linalg as LA
from scipy.spatial.transform import Rotation as R
from specfabpy import specfab as sf
from specfabpy import common as sfcom
from specfabpy import constants as sfconst

L = 8
lm, nlm_len = sf.init(L)
grain_params = sfconst.ice['viscoplastic']['linear']
A, n = 1, 3

x,y,z = np.eye(3)

def evaluate(a2, expname):

    nlm = np.zeros((nlm_len), dtype=np.complex64)
    nlm[:sf.L2len] = sf.a2_to_nlm(a2)
    mi, _ = sfcom.eigenframe(nlm) 
    Eij = sf.Eij_tranisotropic(nlm, *mi, *grain_params) # 3x3 enhancement-factor tensor

    # stress state
    e1, e2 = x, y # coord axes
#    e1, e2 = m2, m3 # eigenframe
    if 0:  tau0 = np.eye(3)/3 - np.einsum('i,j->ij', e1, e1) # compression along e1
    else:  tau0 = np.einsum('i,j->ij', e1, e2) + np.einsum('i,j->ij', e2, e1) # e1-e2 shear
    eps = sf.rheo_fwd_orthotropic(tau0,A,n, *mi, Eij)
    tau = sf.rheo_rev_orthotropic(eps,A,n, *mi, Eij)

    print('\n=== %s ==='%(expname))
    print('\nmi:\n',mi)
    print('\ntau0:\n', tau0)
    print('\neps(tau0):\n', eps)
    print('\ntau(eps(tau0)):\n', tau)
    
    _, ei_tau0 = LA.eig(tau0)
    _, ei_eps = LA.eig(eps)
    
    print('\nei_tau0:\n', ei_tau0)
    print('\nei_eps:\n', ei_eps)
    print('\n*** if above eigenframes ei_* are not equal, then tau and eps are not coaxial ***')
    
    coax = np.einsum('ij,ji',tau0,eps)/(LA.norm(tau0)*LA.norm(eps))
    print('tau0:eps/norm = %.4f'%coax)
    
    
a2_gd = np.diag([0.0,0.5,0.5])
M = R.from_rotvec(np.deg2rad(20) * np.array([0, 1, 0])).as_matrix()
a2 = np.matmul(M.T, np.matmul(a2_gd, M))
evaluate(a2, 'girdle')

a2_sm = np.diag([0.0,0.0,1])
M = R.from_rotvec(np.deg2rad(20) * np.array([0, 1, 0])).as_matrix()
a2 = np.matmul(M.T, np.matmul(a2_sm, M))
evaluate(a2, 'smax')
