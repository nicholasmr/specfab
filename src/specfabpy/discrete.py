#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2023-

"""
Discrete methods
"""

import numpy as np
import scipy.special as sp

from .specfabpy import specfabpy as sf__ # sf private copy 
lm__, nlm_len__ = sf__.init(20)

from . import common as sfcom

def lat2colat(lat, deg=False):   return 90 - lat   if deg else np.pi/2 - lat
def colat2lat(colat, deg=False): return 90 - colat if deg else np.pi/2 - colat

def L2nlmlen(L): return int((L+1)*(L+2)/2)


def cart2sph(v, deg=False):

    """
    Cartesian vector(s) --> spherical coordinates
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


def sphbasis(colat, lon, deg=False):

    """
    Spherical coordinates --> spherical basis vectors
    """

    if deg: colat, lon = np.deg2rad(colat), np.deg2rad(lon)

    rhat = np.array([np.cos(lon)*np.sin(colat), np.sin(lon)*np.sin(colat), +np.cos(colat)]) 
    that = np.array([np.cos(lon)*np.cos(colat), np.sin(lon)*np.cos(colat), -np.sin(colat)]) 
    phat = np.array([-np.sin(lon), np.cos(lon), 0*lon]) 
    
    return (rhat,that,phat) # unit vectors (r, theta, phi)


def idealstressstate(latres=40, lonres=80):

    """
    Idealized stress states w.r.t. spherical unit vectors over all S^2
    
    r = radial vector, t = theta vector, p = phi vector
    """

    vlat = np.deg2rad(np.linspace(-90, 90, latres))
    vlon = np.deg2rad(np.linspace(0, 360, lonres))
    vcolat = lat2colat(vlat, deg=False)

    ### Determine spherical basis vectors
    mlon, mlat = np.meshgrid(vlon, vlat)
    mcolat = lat2colat(mlat, deg=False)
    vr, vt, vp = sphbasis(mcolat, mlon)

    rr = np.einsum('ikl,jkl->ijkl', vr, vr)
    rt = np.einsum('ikl,jkl->ijkl', vr, vt)
    rp = np.einsum('ikl,jkl->ijkl', vr, vp)

    id4 = np.zeros((3,3,latres,lonres)) # repeated identity matrix
    id4[0,0,:,:] = id4[1,1,:,:] = id4[2,2,:,:] = 1
    tau_rr = id4-3*rr
    tau_rt = rt + np.einsum('ijkl->jikl',rt)
    tau_rp = rp + np.einsum('ijkl->jikl',rp)

    return (vr,vt,vp, rr,rt,rp, tau_rr,tau_rt,tau_rp, vlon,vlat, mlon,mlat)
    

"""
Deprecated, use fortran version sf.ri_to_nlm(ri, wi, L)
"""
#def vi2nlm(vi, L=6):

#    """
#    Discrete array of crystallographic axes --> spectral CPO representation
#    """

#    nlm_len = L2nlmlen(L)
#    nlm = np.zeros((nlm_len), dtype=np.complex128) 
#    
#    vi = np.array(vi)
#    N = vi.shape[0] # vi(N,xyz)
#    
#    if L==2:
#        a2 = np.array([ np.einsum('i,j',vi[n],vi[n]) for n in range(N)]).mean(axis=0)
#        nlm[:sf__.L2len] = sf__.a2_to_nlm(a2)    
#    elif L==4:
#        a4 = np.array([ np.einsum('i,j,k,l',vi[n],vi[n],vi[n],vi[n]) for n in range(N)]).mean(axis=0)
#        nlm[:sf__.L4len] = sf__.a4_to_nlm(a4)
#    elif L==6:
#        a6 = np.array([ np.einsum('i,j,k,l,m,n',vi[n],vi[n],vi[n],vi[n],vi[n],vi[n]) for n in range(N)]).mean(axis=0)
#        nlm[:sf__.L6len] = sf__.a6_to_nlm(a6)
#    else:
#        raise ValueError('sfdsc.vi2nlm(): only L <= 6 is supported')

#    return nlm
    
    
def Sl_delta(L):

    nlm_len = L2nlmlen(L)
    nlm = np.zeros((nlm_len), dtype=np.complex128)
    Lrange = np.arange(0,L+1,2) # 0, 2, 4, ...
        
    for ii, (l,m) in enumerate(lm__[:,:nlm_len].T): 
        nlm[ii] = sp.sph_harm(m,l, 0,0) # delta function values
        
    Sl = np.array([sf__.Sl(nlm, l) for l in Lrange])
    Sl /= Sl[0] # normalize
    
    return Sl, Lrange, nlm
    
    
def nlm_ideal_cases(norm=1/np.sqrt(4*np.pi)):

    m = [0,0,1]
    Il24 = [sf__.I20, sf__.I40] # l=2,4, m=0 coefs
    L = 8
    n20_unidir, n40_unidir = np.real(sf__.nlm_ideal(m, 0, L))[Il24]/norm       # delta distributed
    n20_planar, n40_planar = np.real(sf__.nlm_ideal(m, np.pi/2, L))[Il24]/norm # x--y planar distributed
    n20_circle, n40_circle = np.real(sf__.nlm_ideal(m, np.pi/4, L))[Il24]/norm # 45 deg circle distributed

    return ((n20_unidir,n40_unidir), (n20_planar,n40_planar), (n20_circle,n40_circle))
    
    
class DFSE():
    
    """
    Discrete Finite Strain Ellipse (FSE) solver
    
    Based on Lev and Hager (2008) 
    """

    def __init__(self, dim=2, F0=None):
        self.dim = dim
        self.F = F0 if F0 is not None else np.eye(dim)

    def evolve(self, L, dt):
        F0 = self.F.copy()
        A = np.eye(self.dim) - dt/2*L
        B = np.eye(self.dim) + dt/2*L
        self.F = np.matmul(np.linalg.inv(A), np.matmul(B, F0))
        
    def eigenframe(self):
        return sfcom.eigenframe(sfcom.F2C(self.F), modelplane=None) # first entry is largest stretching ratio pair  
   
    
class DDM():

    """
    Discrete Directors Method (DDM) for initially isotropic vector bundle of crystallographic axes
    
    Note: Also implemented as a much faster version in Fortran module, ri_LROT() (dynamics.f90)
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
        self.v += dt*self.velocityfield(self.v, L, iota=self.iota)
        # re-normalize
        norms = np.linalg.norm(self.v, axis=1)
        self.v = np.array([self.v[ii,:]/norms[ii] for ii in np.arange(self.N)])

    def a2(self):
        return np.mean(np.einsum('ni,nj->nij', self.v, self.v), axis=0)

    def eigenframe(self, modelplane=None):
        return sfcom.eigenframe(self.a2(), modelplane=modelplane)

