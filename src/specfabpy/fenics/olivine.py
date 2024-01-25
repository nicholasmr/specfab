#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022-2024

import copy, sys, time, code # code.interact(local=locals())
 
import numpy as np
from dolfin import *

from ..specfabpy import specfabpy as sf__
from .. import constants as sfconst
from .CPO import CPO as fenics_CPO
from .FSE import FSE as fenics_FSE

class OlivineFabric():

    """
    Supports both SDM and FSE representations of CPO (see enable_* flags)
    """

    rho_olivine      = sfconst.olivine['density']
    cij_Abramson1997 = sfconst.olivine['elastic']['Abramson1997']
    cij_Jacobsen2008 = sfconst.olivine['elastic']['Jacobsen2008']

    def __init__(self, meshargs, L=8, fabrictype='A', nu_realspace=1e-3, nu_multiplier=1, enable_SDM=True, enable_FSE=False, Cij=None, rho=None): 
        mesh, boundaries, ds, n = meshargs
        self.L = L
        self.fabrictype = fabrictype
        sf__.init(self.L)
        self.set_elastic_params(Cij=Cij, rho=rho)
        self.enable_SDM = enable_SDM
        self.enable_FSE = enable_FSE
        if self.enable_FSE: self.FSE = fenics_FSE(mesh, boundaries, ds, n)
                
        kwargs = dict(modelplane='xz', nu_multiplier=nu_multiplier, nu_realspace=nu_realspace)
        self.SDM_n = fenics_CPO(mesh, boundaries, ds, n, self.L, **kwargs)
        self.SDM_b = fenics_CPO(mesh, boundaries, ds, n, self.L, **kwargs)
        self.SDM_n.initialize(wr=None) # isotropic
        self.SDM_b.initialize(wr=None) 
        self.SDM_n.set_BCs([], [], []) # no BCs
        self.SDM_b.set_BCs([], [], [])
        
        if fabrictype not in ['A', 'B', 'C']:
            raise ValueError('Fabric type "%s" not among supported types "A", "B", "C"'%(fabrictype))

    def set_state(self, sb, sn, interp=True):
        if interp:
            raise ValueError('OlivineFabric.set_state() supports only setting function space vars, not interpolating expressions or constants.')
        else:
            self.SDM_b.w.assign(sb)
            self.SDM_n.w.assign(sn)
            self.SDM_b.w_prev.assign(sb)
            self.SDM_n.w_prev.assign(sn)

    def get_state(self, x, y):
        return (self.SDM_n.get_nlm(x,y), self.SDM_b.get_nlm(x,y))

    def get_FSE(self, x, y):
        return self.FSE.eigenframe(x,y)
        
    def evolve(self, u, dt, iota_n=+1, iota_b=-1, Gamma0=None, Lambda0=None):
        if self.enable_SDM:
            self.SDM_n.evolve(u, dt, iota=iota_n, Gamma0=None, Lambda0=None)
            self.SDM_b.evolve(u, dt, iota=iota_b, Gamma0=None, Lambda0=None)
        if self.enable_FSE: 
            self.FSE.evolve(u, dt)
            
    def set_elastic_params(self, Cij=None, rho=None):
        if Cij is not None: 
            # Lame_grain = (lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3) w.r.t. mi' axes
            l = np.concatenate(([None], sf__.Cij_to_Lame_orthotropic(Cij))) # prepend "None" for easy indexing
            # (1,2,3, 4,5,6, 7,8,9)
            # Lame_grain should be (lam_bb,lam_nn,lam_vv, lam_nv,lam_bv,lam_nb, mu_b,mu_n,mu_v)
            if self.fabrictype == 'A': l = l[1:] # already correct order since (b,n,v)=(m1',m2',m3') for A-type
            if self.fabrictype == 'B': l = sf.Lame_olivine_A2X(l,'B') #[l[3],l[2],l[1], l[6],l[5],l[4], l[9],l[8],l[7]] # (b,n,v)=(m3',m2',m1')
            if self.fabrictype == 'C': l = sf.Lame_olivine_A2X(l,'C') #[l[3],l[1],l[2], l[6],l[4],l[5], l[9],l[7],l[8]] # (b,n,v)=(m3',m1',m2')
            self.Lame_grain = l
        if rho is not None: 
            self.rho = rho

    def get_elastic_velocities(self, x,y, theta,phi, alpha=1):
        alpha = 1 # only strain homogenization is so far supported (alpha=1)
        nlm, blm = self.get_state(x,y)
        vlm = 0*nlm # estimate from joint ODF (nlm and blm) if zero
        vS1, vS2, vP = sf__.Vi_elastic_orthotropic(blm, nlm, vlm, alpha,self.Lame_grain,self.rho, theta,phi) # calculate elastic phase velocities using specfab
        return (vP, vS1, vS2)

