#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022-2024

import copy, sys, time, code # code.interact(local=locals())
 
import numpy as np
from dolfin import *

from ..specfabpy import specfabpy as sf__
from .. import constants as sfconst
from .CPO import CPO as fenics_CPO
from .FSE import FSE as fenics_FSE
from .enhancementfactor import *

class OlivineFabric():

    """
    Supports both SDM and FSE representations of CPO (see enable_* flags)
    """

    rho_olivine      = sfconst.olivine['density']
    cij_Abramson1997 = sfconst.olivine['elastic']['Abramson1997']
    cij_Jacobsen2008 = sfconst.olivine['elastic']['Jacobsen2008']

    def __init__(self, mesh, boundaries, L=8, nu_realspace=1e-3, nu_multiplier=1, modelplane='xz', ds=None, nvec=None, fabrictype='A', enable_SDM=True, enable_FSE=False, \
                    Eij_grain=(1,1,1, 1,1,1), alpha=0, n_grain=1, \
                    Cij=None, rho=None): 
                    
        self.mesh, self.boundaries = mesh, boundaries
        self.L = L
        self.enable_SDM = enable_SDM
        self.enable_FSE = enable_FSE
        self.fabrictype = fabrictype # only spatio-temporal constant type allowed for now
        
        if self.fabrictype not in ['A', 'B', 'C']:
            raise ValueError('Fabric type "%s" not among supported types "A", "B", "C"'%(fabrictype))

        for key in ['SDM_n', 'SDM_b']:
            SDM = fenics_CPO(self.mesh, self.boundaries, self.L, modelplane=modelplane, nu_multiplier=nu_multiplier, nu_realspace=nu_realspace)
            SDM.initialize(sr=None) # isotropic
            SDM.set_BCs([], [], []) # no BCs
            setattr(self, key, SDM)
        
        if self.enable_FSE: 
            self.FSE = fenics_FSE(mesh, boundaries)

        self.sf = sf__
        self.lm, self.nlm_len = self.sf.init(self.L)
        self.set_elastic_params(Cij=Cij, rho=rho)

        if alpha > 0:    raise ValueError('Only Sachs (alpha = 1) homogenization is supported.')
        if n_grain != 1: raise ValueError("Only linear-viscous grain rheology (n'=1) is supported.")
        # Grain enhancement factors w.r.t. {q1,q2,q3} = {b,n,v} axes.
        # Set Enb >> 1 to make n--b slip system comparatively soft.
        # E.g.: self.Eij_grain = (1,1,1, 1,1,Enb) # Voigt order: (Ebb, Enn, Evv, Env, Ebv, Enb)
        self.grain_params = (Eij_grain, alpha, n_grain)
        self.enhancementfactor = EnhancementFactor(self.mesh, self.L, modelplane=modelplane)
        self.update_Eij()
                

    def set_state(self, sb, sn, interp=True):
        if interp:
            raise ValueError('OlivineFabric.set_state() supports only setting function space vars, not interpolating expressions or constants.')
        else:
            self.SDM_b.s.assign(sb)
            self.SDM_n.s.assign(sn)
            self.SDM_b.s0.assign(sb)
            self.SDM_n.s0.assign(sn)

    def get_state(self, x, y):
        return (self.SDM_n.get_nlm(x,y), self.SDM_b.get_nlm(x,y))

    def get_FSE(self, x, y):
        return self.FSE.eigenframe(x,y)
        
    def evolve(self, u, S, dt, iota_n=+1, iota_b=-1, Gamma0=None, Lambda0=None):
        if self.enable_SDM:
            self.SDM_n.evolve(u, S, dt, iota=iota_n, Gamma0=Gamma0, Lambda0=Lambda0)
            self.SDM_b.evolve(u, S, dt, iota=iota_b, Gamma0=Gamma0, Lambda0=Lambda0)
            self.update_Eij()
        if self.enable_FSE: 
            self.FSE.evolve(u, dt)
            
    def update_Eij(self):
        self.mi, self.Eij, self.ai = self.enhancementfactor.Eij_orthotropic(self.SDM_b.s, self.SDM_n.s, *self.grain_params, ei=())
        self.xi, self.Exij, _      = self.enhancementfactor.Eij_orthotropic(self.SDM_b.s, self.SDM_n.s, *self.grain_params, ei=np.eye(3)) 
        # ... unpack
        self.m1, self.m2, self.m3 = self.mi # eigenvectors (presumed fabric and rheological symmetry directions)
        self.E11, self.E22, self.E33, self.E23, self.E31, self.E12 = self.Eij  # eigenenhancements
        self.Exx, self.Eyy, self.Ezz, self.Eyz, self.Exz, self.Exy = self.Exij # Cartesian enhancements
        self.a1, self.a2, self.a3 = self.ai # CPO eigenvalues
            
    def set_elastic_params(self, Cij=None, rho=None):
        if Cij is not None: 
            # Lame_grain = (lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3) w.r.t. mi' axes
            l = np.concatenate(([None], self.sf.Cij_to_Lame_orthotropic(Cij))) # prepend "None" for easy indexing
            # (1,2,3, 4,5,6, 7,8,9)
            # Lame_grain should be (lam_bb,lam_nn,lam_vv, lam_nv,lam_bv,lam_nb, mu_b,mu_n,mu_v)
            if self.fabrictype == 'A': l = l[1:] # already correct order since (b,n,v)=(m1',m2',m3') for A-type
            if self.fabrictype == 'B': l = self.sf.Lame_olivine_A2X(l,'B') #[l[3],l[2],l[1], l[6],l[5],l[4], l[9],l[8],l[7]] # (b,n,v)=(m3',m2',m1')
            if self.fabrictype == 'C': l = self.sf.Lame_olivine_A2X(l,'C') #[l[3],l[1],l[2], l[6],l[4],l[5], l[9],l[7],l[8]] # (b,n,v)=(m3',m1',m2')
            self.Lame_grain = l
        if rho is not None: 
            self.rho = rho
            
    def get_elastic_velocities(self, x,y, theta,phi, alpha=1):
        alpha = 1 # only strain homogenization is so far supported (alpha=1)
        nlm, blm = self.get_state(x,y)
        vlm = 0*nlm # estimate from joint ODF (nlm and blm) if zero
        vS1, vS2, vP = self.sf.Vi_elastic_orthotropic(blm, nlm, vlm, alpha,self.Lame_grain,self.rho, theta,phi) # calculate elastic phase velocities using specfab
        return (vP, vS1, vS2)
        
    def get_elastic_velocities__isotropic(self, alpha=1):
        alpha = 1
        theta,phi = 0,0
        nlm = blm = vlm = [1/np.sqrt(4*np.pi)] + [0]*(self.nlm_len-1)
        vS1, vS2, vP = self.sf.Vi_elastic_orthotropic(blm, nlm, vlm, alpha,self.Lame_grain,self.rho, theta,phi)
        return (vP, vS1, vS2)
    

