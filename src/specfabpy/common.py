#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2024

"""
Library of commonly used routines
"""

import numpy as np
import copy, sys, os, code # code.interact(local=locals())
from .specfabpy import specfabpy as sf__ # sf private copy 
lm__, nlm_len__ = sf__.init(20)

def eigenframe(Z, modelplane=None, symframe=-1):

    """
    Principal fabric frame (typically a^(2) principal axes).
    
    if Z[nn,3,3] --> list of a2 matrices
    if Z[3,3]    --> single  a2 matrix
    if Z[nn,:]   --> list of nlm vectors
    if Z[:]      --> single  nlm vector    
    """
    
    if Z.shape[-1] > 3: 
        # nlm vector passed, or list thereof was passed
        nlm = np.array(Z, ndmin=2) # ensure nlm[nn,coef] format even if nlm[coef] was passed 
        Q = np.array([sf__.a2(nlm[nn,:]) if (symframe == -1) else sf__.a4_eigentensors(nlm[nn,:])[symframe] for nn in range(nlm.shape[0])])
    else:     
        # assume Z = a^(2) passed directly, or list thereof
        Z = np.array(Z, ndmin=3) # ensure Z[nn,3,3] format even if Z[3,3] matrix was passed
        Q = Z

    if modelplane is None: modelplane='ij'
    ei, lami = sf__.eigframe_arr(Q, modelplane)
#    ei = np.moveaxis(ei, -1, -2) # ei[node,i,xyz] --> ei[node,xyz,i]

    return (ei[0], lami[0]) if Q.shape[0]==1 else (ei, lami) # ei[node,i,xyz], lami[node,i]

    
def xi_tile(N):     return ei_tile(np.eye(3), N)
def ei_tile(ei, N): return (np.tile(ei[0],(N,1)), np.tile(ei[1],(N,1)), np.tile(ei[2],(N,1)))
    
    
def pfJ(nlm): 

    """
    Pole figure J (pfJ) index
    """

    k = np.sqrt(4*np.pi)**2
    if len(nlm.shape) == 2: J_i = [ k*np.sum(np.square(np.absolute(nlm[ii]))) for ii in range(nlm.shape[0]) ]
    else:                   J_i =   k*np.sum(np.square(np.absolute(nlm)))
    return J_i
    
   
def F2C(F):

    """
    Deformation gradient to C
    """
    
    Finv = np.linalg.inv(F)
    C = np.dot(Finv.T, Finv)
    #C = np.dot(F.T, F) # wrong principal directions?
    return C
    
    
def mat2d(D3, modelplane):

    """
    Model-plane strain-rate tensor from R^3 strain-rate tensor
    """
    
    if modelplane=='xy':
        D2 = np.array([ \
            [D3[0,0], D3[0,1]], \
            [D3[1,0], D3[1,1]] \
        ])
    
    elif modelplane=='xz':  
        D2 = np.array([ \
            [D3[0,0], D3[0,2]], \
            [D3[2,0], D3[2,2]] \
        ])
        
    else:
        raise ValueError('invalid modelplane "%s"'%(modelplane))
        
    return D2
    
    
def mat3d_arr(D2, modelplane, reshape=False):

    """
    Overloaded version of mat3d() that can take arrays of D2
    """   
     
    # 2D to 3D strain-rate/stress tensor assuming (i) tr=0 and (ii) out-of-modelplane shear components vanish
    if len(D2.shape) == 1: D2 = np.array(D2, ndmin=2) # generalize so that this routine works for arrays of matrices, too
    D3 = np.array([ mat3d(D2[nn], modelplane, reshape=reshape) for nn in range(D2.shape[0]) ])
    return D3[0] if D3.shape[0]==1 else D3
    
    
def mat3d(D2, modelplane, reshape=False): 

    """
    R^3 strain-rate tensor from model-plane strain-rate tensor
    """

    if reshape: D2 = np.reshape(D2,(2,2)) # reshape if vectorized
    trace = D2[0,0] + D2[1,1]
    
    if modelplane=='xy':
        D3 = np.array([ [D2[0,0], D2[0,1], 0],\
                        [D2[1,0], D2[1,1], 0],\
                        [0,             0,  -trace] ] )  
    
    elif modelplane=='xz':  
        D3 = np.array([ [D2[0,0],      0,    D2[0,1]],\
                        [0,       -trace,    0      ],\
                        [D2[1,0],      0,    D2[1,1]] ] )  
                        
    else:
        raise ValueError('invalid modelplane "%s"'%(modelplane))
        
    return D3
    
