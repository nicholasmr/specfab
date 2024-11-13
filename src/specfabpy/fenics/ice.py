#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>

import copy, sys, time, code # code.interact(local=locals())
import numpy as np
import matplotlib.tri as tri

from .. import constants as sfconst
from .CPO import *
from .enhancementfactor import *

class IceFabric(CPO):

    def __init__(self, mesh, boundaries, *args, L=8, nu_realspace=1e-3, nu_multiplier=1, modelplane='xz', ds=None, nvec=None, \
                        Eij_grain=(1,1), alpha=0, n_grain=1, symframe=-1, \
                        Cij=sfconst.ice['elastic']['Bennett1968'], rho=sfconst.ice['density'], **kwargs): 

        super().__init__(mesh, boundaries, L, nu_multiplier=nu_multiplier, nu_realspace=nu_realspace, modelplane=modelplane, ds=ds, nvec=nvec)
        self.initialize(sr=None) # isotropic
        self.set_BCs([], [], []) # no BCs

        self.grain_params = (Eij_grain, alpha, n_grain)
        self.Lame_grain = self.sf.Cij_to_Lame_tranisotropic(Cij) 
        self.rho = rho

        self.symframe = symframe
        self.enhancementfactor = EnhancementFactor(mesh, L, symframe=symframe, modelplane=modelplane)
        self.update_Eij()
                
    def get_state(self, *args, **kwargs): 
        return self.get_nlm(*args, **kwargs) # alias
        
    def set_state(self, *args, **kwargs):
        super().set_state(*args, **kwargs)
        self.update_Eij()
        
    def evolve(self, *args, **kwargs):
        super().evolve(*args, **kwargs)
        self.update_Eij()

    def solvesteady(self, u, S, **kwargs):
        super().solvesteady(u, S, **kwargs)
        self.update_Eij()
        
    def update_Eij(self):
        self.mi, self.Eij, self.lami = self.enhancementfactor.Eij_tranisotropic(self.s, *self.grain_params, ei=())
        self.xi, self.Exij, _        = self.enhancementfactor.Eij_tranisotropic(self.s, *self.grain_params, ei=np.eye(3))
        self.pfJ = self.enhancementfactor.pfJ(self.s)
        # ... unpack
        self.m1, self.m2, self.m3 = self.mi # rheological symmetry directions
        self.E11, self.E22, self.E33, self.E23, self.E31, self.E12 = self.Eij  # eigenenhancements
        self.Exx, self.Eyy, self.Ezz, self.Eyz, self.Exz, self.Exy = self.Exij # Cartesian enhancements
        self.lam1, self.lam2, self.lam3 = self.lami # fabric eigenvalues \lambda_i (not a2 eigenvalues unless symframe=-1)
            
    def E_CAFFE(self, u, *args, **kwargs):
        return self.enhancementfactor.E_CAFFE(self.s, u, *args, **kwargs)
            
    def E_EIE(self, u, *args, **kwargs):
        return self.enhancementfactor.E_EIE(u, self.Eij, self.mi, *args, **kwargs)
            
    def get_elastic_velocities(self, x,y, theta,phi, alpha=1):
        nlm = self.get_state(x,y)
        vS1, vS2, vP = sf__.Vi_elastic_tranisotropic(nlm, alpha,self.Lame_grain,self.rho, theta,phi) # calculate elastic phase velocities using specfab
        return (vP, vS1, vS2)
        
    def Gamma0_Lilien(self, u, T, A=4.3e7, Q=3.36e4):
        # DDRX rate factor from Dome C ice-core calibration experiment (Lilien et al., 2023, p. 7)
        R = Constant(8.314) # gas constant (J/mol*K)
        D = sym(grad(u))
        epsE = sqrt(inner(D,D)/2)
        return project(epsE*Constant(A)*exp(-Constant(Q)/(R*T)), self.R)
        
    def df2np(self, F, withcoords=False):
        coords = copy.deepcopy(self.mesh.coordinates().reshape((-1, 2)).T)
        triang = tri.Triangulation(*coords, triangles=self.mesh.cells())    
        Q = FunctionSpace(self.mesh, 'CG', 2)
        F_np = project(F,Q).compute_vertex_values(self.mesh)
        return (triang, F_np, coords) if withcoords else (triang, F_np)
        
