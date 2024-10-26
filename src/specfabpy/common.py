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
    
    if len(Z.shape) >= 1:
        # if Z is multidimensional then assume nlm[nn,coefs] was passed
        nlm = np.array(Z, ndmin=2) # ensure nlm[nn,coef] format even if nlm[coef] was passed 
        Q = np.array([sf__.a2(nlm[nn,:]) if (symframe == -1) else sf__.a4_eigentensors(nlm[nn,:])[symframe] for nn in range(nlm.shape[0])])
    else:                 
        Z = np.array(Z, ndmin=2)
        Q = Z # assume Z = a^(2) passed directly
        
    lami, ei = np.linalg.eig(Q)

    if modelplane is None: 
        I = np.flip(lami.argsort()) # largest eigenvalue pair is first entry
        ei, lami = ei[:,:,I], lami[:,I]
    else:
        kk = 1 if modelplane=='xz' else 2 # xz else xy
        for nn in range(lami.shape[0]):
            lami_, ei_ = lami[nn,:], ei[nn,:,:] # node nn
            I = abs(ei_[kk,:]).argsort() # last entry is ei with largest y-value (component out of xz modelplane), or z value for xy modelplane
            II = np.flip(lami_[I[0:2]].argsort()) # largest eigenvalue sorting of *in-plane ei vectors*
            I[0:2] = I[II] # rearrange in-plane eigenpairs so that largest eigenvalue pair is first
            lami[nn,:], ei[nn,:,:] = lami_[I], ei_[:,I]
        
    return (ei[0], lami[0]) if Z.shape[0]==1 else (ei, lami) # [node], [, xyz comp], eigenpair
    
    
def pfJ(nlm): 

    """
    Pole figure J (pfJ) index
    """
    # @TODO: nlm index ordering should be transposed to match [node, comp] ordering of other routines
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
    
