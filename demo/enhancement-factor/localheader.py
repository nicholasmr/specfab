#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, code # code.interact(local=locals())

import numpy as np
from specfabpy import specfab as sf__

# private
L__ = 8
lm__, nlm_len__ = sf__.init(L__) 

m, t = np.array([0,0,1]), np.array([1,0,0])
p, q = (m+t)/np.sqrt(2), (m-t)/np.sqrt(2)
mm, mt, pq = np.tensordot(m,m, axes=0), np.tensordot(m,t, axes=0), np.tensordot(p,q, axes=0)

tau0 = 1 # results are independant of the stress magnitude
tau_mm = tau0*(np.identity(3)/3 - mm) 
tau_mt = tau0*(mt + mt.T) 
tau_pq = tau0*(pq + pq.T)

def Eij_maps_tranisotropic(X,Y, paramcombo, Ecc=None, Eca=None, alpha=None, n_grain=1):

    X_map, Y_map = np.meshgrid(X, Y, indexing='xy')

    size = (len(Y),len(X))
    Emm_map, Emt_map, Epq_map = np.zeros(size), np.zeros(size), np.zeros(size)

    nlm_unidir = sf__.nlm_ideal(m, 0, L__)

    for ii,y in enumerate(Y):
        for jj,x in enumerate(X):
        
            if paramcombo == 'Taylor':
                Ecc   = x
                Eca   = y                
                alpha = 1
            elif paramcombo == 'Sachs':
                Ecc   = x
                Eca   = y            
                alpha = 0
            elif paramcombo == 'Mixed':
                alpha = x
                Eca   = y
            else:
                raise ValueError('argument "paramcombo" not understood')
                
            grain_params = ((Ecc,Eca), alpha, n_grain)
            
            Emm_map[ii,jj] = sf__.Evw_tranisotropic(nlm_unidir, m,m,tau_mm, *grain_params)
            Emt_map[ii,jj] = sf__.Evw_tranisotropic(nlm_unidir, m,t,tau_mt, *grain_params)
            Epq_map[ii,jj] = sf__.Evw_tranisotropic(nlm_unidir, p,q,tau_pq, *grain_params)
            
            
    return Emm_map, Emt_map, Epq_map, X_map, Y_map
    
    
def Eij_maps_orthotropic(X,Y, paramcombo, Eii_grain=None, n_grain=1):

    X_map, Y_map = np.meshgrid(X, Y, indexing='xy')

    size = (len(Y),len(X))
    Emm_map, Emt_map, Epq_map = np.zeros(size), np.zeros(size), np.zeros(size)

    Eij_grain = np.ones(6)

    nlm_2 = sf__.nlm_ideal(m, 0, L__) # n unidir in m dir.
    nlm_1 = sf__.nlm_ideal(t, 0, L__) # b unidir in t dir.
    nlm_3 = 0*nlm_1 # infer from nlm_1 and nlm_2

    for ii,y in enumerate(Y):
        for jj,x in enumerate(X):
        
            if paramcombo == 'Taylor':
                Eij_grain[[0,1,2]] = x
                Eij_grain[5] = y # x-y shear
                alpha = 1
            elif paramcombo == 'Sachs':
                Eij_grain[[0,1,2]] = x
                Eij_grain[5] = y # x-y shear
                alpha = 0
            elif paramcombo == 'Mixed':
                Eij_grain[[0,1,2]] = Eii_grain
                Eij_grain[5] = y # x-y shear
                alpha = x
            else:
                raise ValueError('argument "paramcombo" not understood')
                
            grain_params = (Eij_grain, alpha, n_grain)
            
            Emm_map[ii,jj] = sf__.Evw_orthotropic(nlm_1,nlm_2,nlm_3, m,m,tau_mm, *grain_params)
            Emt_map[ii,jj] = sf__.Evw_orthotropic(nlm_1,nlm_2,nlm_3, m,t,tau_mt, *grain_params)
            Epq_map[ii,jj] = sf__.Evw_orthotropic(nlm_1,nlm_2,nlm_3, p,q,tau_pq, *grain_params)
            
    return Emm_map, Emt_map, Epq_map, X_map, Y_map

