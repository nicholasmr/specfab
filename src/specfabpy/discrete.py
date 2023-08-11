#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2023

"""
Discrete methods
"""

import numpy as np
import scipy.special as sp


def lat2colat(lat, deg=False):   return 90 - lat   if deg else np.pi/2 - lat
def colat2lat(colat, deg=False): return 90 - colat if deg else np.pi/2 - colat

def L2nlmlen(L): return int((L+1)*(L+2)/2)

def cart2sph(v, deg=False):

    """
    Cartesian vector(s) <-- spherical coordinates
    """

    v = np.array(v)
    if len(v.shape) == 1: x,y,z = v
    else:                 x,y,z = v[:,0],v[:,1],v[:,2]

    lon = np.arctan2(y,x)
    colat = np.arctan2(np.sqrt(np.power(x,2) + np.power(y,2)), z)
    lat = np.pi/2 - colat
    return (np.rad2deg(lat), np.rad2deg(colat), np.rad2deg(lon)) if deg else (lat, colat, lon)
    

def sph2cart(colat, lon, deg=False):

    """
    Spherical coordinates --> Cartesian vector(s)
    """
    
    if deg: colat, lon = np.deg2rad(colat), np.deg2rad(lon)
    colat, lon = np.atleast_1d(colat), np.atleast_1d(lon)
    
    r = np.array([ [np.cos(lon[ii])*np.sin(colat[ii]), np.sin(lon[ii])*np.sin(colat[ii]), +np.cos(colat[ii])] for ii in range(len(colat)) ]) 
    return r


def vi2nlm(vi, sf, L=6):

    """
    Discrete array of crystallographic axes --> spectral CPO representation
    """

    nlm_len = L2nlmlen(L)
    nlm = np.zeros((nlm_len), dtype=np.complex128) 

    if L==2:
        a2 = np.array([ np.einsum('i,j',v,v) for v in vi]).mean(axis=0)
        nlm[:sf.L2len] = sf.a2_to_nlm(a2)    
    elif L==4:
        a4 = np.array([ np.einsum('i,j,k,l',v,v,v,v) for v in vi]).mean(axis=0)
        nlm[:sf.L4len] = sf.a4_to_nlm(a4)
    elif L==6:
        a6 = np.array([ np.einsum('i,j,k,l,m,n',v,v,v,v,v,v) for v in vi]).mean(axis=0)
        nlm[:sf.L6len] = sf.a6_to_nlm(a6)
    else:
        raise ValueError('sfdsc.vi2nlm(): only L <= 6 is supported')

    return nlm


def sphericalbasisvectors(colat, lon, deg=False):

    """
    Spherical coordinates --> spherical basis vectors
    """

    if deg: colat, lon = np.deg2rad(colat), np.deg2rad(lon)

    rhat = np.array([np.cos(lon)*np.sin(colat), np.sin(lon)*np.sin(colat), +np.cos(colat)]) 
    that = np.array([np.cos(lon)*np.cos(colat), np.sin(lon)*np.cos(colat), -np.sin(colat)]) 
    phat = np.array([-np.sin(lon), np.cos(lon), 0*lon]) 
    
    return (rhat,that,phat) # unit vectors (r, theta, phi)


def Sl_delta(lm, sf):

    L = lm[0,-1]
    nlm_len = L2nlmlen(L)
    nlm = np.zeros((nlm_len), dtype=np.complex128)
    Lrange = np.arange(0,L+1,2) # 0, 2, 4, ...
        
    for ii, (l,m) in enumerate(lm.T): 
        nlm[ii] = sp.sph_harm(m,l, 0,0) # delta function values
        
    Sl = np.array([sf.Sl(nlm, l) for l in Lrange])
    Sl /= Sl[0] # normalize
    
    return Sl, Lrange, nlm
    

class DDM():

    """
    Discrete Directors Method (DDM) for initially isotropic vector bundle of crystallographic axes
    """

    def __init__(self, iota=+1, N=1000, v0=None):
    
        self.iota = iota
        self.N = N
        self.v = np.zeros((N,3))        

        # Initialize with randomly distributed vector bundle
        if v0 is None:
            for ii in np.arange(N):
                phi = np.random.uniform(0, 2*np.pi)
                theta = np.arccos(1-2*np.random.uniform(0, 1))
                x = np.sin(theta) * np.cos(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(theta)
                self.v[ii,:] = [x,y,z]
                
        # If more than 3 vectors, set the passed list of init vectors
        elif len(v0.flatten()) > 3:
            self.N = v0.shape[0]
            self.v = v0 # (num, xyz)
            
        else:
            for ii in np.arange(self.N): self.v[ii,:] = v0


    @staticmethod 
    def velocityfield(v0, L, iota=1):
    
        D = (L + L.T)/2 # L is velocity gradient
        W = (L - L.T)/2
        
        v0sq = np.einsum('ni,nj->nij', v0,v0) 
        dvdt = np.einsum('ij,nj->ni', W, v0) # W.v
        asymprod = np.einsum('nij,jk->nik', v0sq, D) - np.einsum('ij,njk->nik', D, v0sq)
        dvdt += iota*np.einsum('nij,nj->ni', asymprod, v0) # iota*(v^2.D - D.v^2).v

        return dvdt


    def evolve(self, L, dt):
        
        self.v += dt*self.velocityfield(self.v, L)
        # re-normalize
        norms = np.linalg.norm(self.v, axis=1)
        self.v = np.array([self.v[ii,:]/norms[ii] for ii in np.arange(self.N)])


    def get_eigenframe(self):
    
        a2_n = np.einsum('ni,nj->nij', self.v, self.v)
        a2 = np.mean(a2_n, axis=0)
        print(a2)
        eigvals, eigvecs = np.linalg.eig(a2)
        I = eigvals.argsort() 
        return eigvecs[:,I], eigvals[I] # last entry is largest eigenvalue pair
        
    
    def get_principal_dir_2D(self, outofplane='y'):

        eigvecs, eigvals = self.get_eigenframe()
#        print('--------')
#        print(eigvals)
#        print(eigvecs)
        if outofplane == 'y': 
            e1 = eigvecs[[0,2],-1] # x,z components of large eigenvalue vector
            e1 /= np.linalg.norm(e1)
            return e1
        else: return None

