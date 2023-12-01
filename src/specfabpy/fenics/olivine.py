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

    def __init__(self, mesh, boundaries, ds, n, L, nu_realspace=1e-3, nu_multiplier=1, enable_SDM=True, enable_FSE=False): 

        self.L = L
        sf__.init(self.L)
    
        self.enable_SDM = enable_SDM
        self.enable_FSE = enable_FSE
        
        kwargs_CPOmodel = dict(modelplane='xz', nu_multiplier=nu_multiplier, nu_realspace=nu_realspace)
        self.SDM_n = fenics_CPO(mesh, boundaries, ds, n, self.L, **kwargs_CPOmodel)
        self.SDM_b = fenics_CPO(mesh, boundaries, ds, n, self.L, **kwargs_CPOmodel)
        self.SDM_n.initialize(wr=None) # isotropic
        self.SDM_b.initialize(wr=None) 
        self.SDM_n.set_BCs([], [], []) # no BCs
        self.SDM_b.set_BCs([], [], [])
        
        if self.enable_FSE: 
            self.FSE = fenics_FSE(mesh, boundaries, ds, n)

    
    def set_elastic_params(self, Cij=None, rho=None):
        if Cij is not None: self.Lame_grain = sf__.Cij_to_Lame_orthotropic(Cij) 
        if rho is not None: self.rho = rho


    def get_statevec(self, x, y):
        return self.SDM_n.get_nlm(x,y), self.SDM_b.get_nlm(x,y)


    def get_FSE(self, x, y):
        return self.FSE.eigenframe(x,y)
        
        
    def evolve(self, u, dt, iota_n=+1, iota_b=-1):
        if self.enable_SDM:
            self.SDM_n.evolve(u, dt, iota=iota_n, Gamma0=None, Lambda0=None)
            self.SDM_b.evolve(u, dt, iota=iota_b, Gamma0=None, Lambda0=None)
        if self.enable_FSE: 
            self.FSE.evolve(u, dt)


    def get_seismicvelocities(self, x,y, theta,phi):
        nlm, blm = self.get_statevec(x,y)
        vlm = 0*nlm # estimate from joint ODF (nlm and blm) if zero
        alpha = 1 # only strain homogenization is so far supported (alpha=1)
        vS1, vS2, vP = sf__.Vi_elastic_orthotropic(blm, nlm, vlm, alpha,self.Lame_grain,self.rho, theta,phi) # calculate elastic phase velocities using specfab
        return (vP, vS1, vS2)

