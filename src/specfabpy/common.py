#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2024

"""
Library of commonly used routines
"""

import numpy as np
from .specfabpy import specfabpy as sf__ # sf private copy 
lm__, nlm_len__ = sf__.init(20)


def eigenframe(Z, modelplane=None, symframe=-1):

    """
    Principal fabric frame (typically a^(2) principal axes).
    """

    # Is input one-dimensional (vector)? => assume nlm was passed
    if len(Z.shape) == 1: 
        Q = sf__.a2(Z) if (symframe == -1) else sf__.a4_eigentensors(Z)[symframe]
    else:                 
        Q = Z # assume Z = a^(2) passed directly
        
    eigvals, eigvecs = np.linalg.eig(Q)
    
    if modelplane is None: 
        I = np.flip(eigvals.argsort()) # largest eigenvalue pair is first entry
        
    else:
        if   modelplane=='xy': I = abs(eigvecs[2,:]).argsort() # last entry is eigenvector with largest z-value (component out of model plane)
        elif modelplane=='xz': I = abs(eigvecs[1,:]).argsort() # ... same but for y-value
        
        if eigvals[I[0]] < eigvals[I[1]]: 
            I[[0,1]] = I[[1,0]] # swap sorted indices so largest eigenvalue entry is first
    
    return eigvecs[:,I], eigvals[I]
    
    
def pfJ(nlm): 

    """
    Pole figure J (pfJ) index
    """
    
    k = np.sqrt(4*np.pi)**2
    if len(nlm.shape) == 2: J_i = [ k*np.sum(np.square(np.absolute(nlm[:,ii]))) for ii in range(nlm.shape[1]) ]
    else:                   J_i =   k*np.sum(np.square(np.absolute(nlm[:,ii])))
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
    
