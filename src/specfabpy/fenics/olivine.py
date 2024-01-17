#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

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

    def __init__(self, meshargs, L=8, nu_realspace=1e-3, nu_multiplier=1, enable_SDM=True, enable_FSE=False, Cij=None, rho=None): 
        mesh, boundaries, ds, n = meshargs
        self.L = L
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
        if Cij is not None: self.Lame_grain = sf__.Cij_to_Lame_orthotropic(Cij) 
        if rho is not None: self.rho = rho

    def get_elastic_velocities(self, x,y, theta,phi, alpha=1):
        alpha = 1 # only strain homogenization is so far supported (alpha=1)
        nlm, blm = self.get_state(x,y)
        vlm = 0*nlm # estimate from joint ODF (nlm and blm) if zero
        vS1, vS2, vP = sf__.Vi_elastic_orthotropic(blm, nlm, vlm, alpha,self.Lame_grain,self.rho, theta,phi) # calculate elastic phase velocities using specfab
        return (vP, vS1, vS2)

