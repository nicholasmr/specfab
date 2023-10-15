#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

"""
Library of commonly used routines
"""

import numpy as np

def eigenframe(a2, modelplane=None):

    """
    Principal frame
    """

    eigvals, eigvecs = np.linalg.eig(a2)
    
    if modelplane is None: 
        I = np.flip(eigvals.argsort()) # largest eigenvalue pair is first entry
        
    else:
        if   modelplane=='xy': I = abs(eigvecs[2,:]).argsort() # last entry is eigenvector with largest z-value (component out of model plane)
        elif modelplane=='xz': I = abs(eigvecs[1,:]).argsort() # ... same but for y-value
        
        if eigvals[I[0]] < eigvals[I[1]]: 
            I[[0,1]] = I[[1,0]] # swap sorted indices so largest eigenvalue entry is first
    
    return eigvecs[:,I], eigvals[I]
    
    
def F2C(F):

    """
    Deformation gradient to C
    """
    
    Finv = np.linalg.inv(F)
    C = np.dot(Finv.T, Finv)
    #C = np.dot(F.T, F) # wrong principal directions?
    return C
    
    
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
