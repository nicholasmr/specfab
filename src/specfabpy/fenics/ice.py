#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>

import copy, code # code.interact(local=locals())
import numpy as np
import matplotlib.tri as tri

from .. import constants as sfconst
from .CPO import *
from .enhancementfactor import *

class IceFabric(CPO):

    def __init__(
        self, mesh, boundaries, L, *args, 
        nu_multiplier=1, nu_realspace=1e-3, modelplane='xz', symframe=-1, 
        enhancementmodel='LTS', Eij_grain=(1,1), alpha=0, n_grain=1, CAFFE_params=(0.1, 10), n_EIE=3, 
        ds=None, nvec=None, setextra=True, 
        Cij=sfconst.ice['elastic']['Bennett1968'], rho=sfconst.ice['density'], **kwargs
    ): 

        ### Check args
        # ...done in CPO() and EnhancementFactor()

        ### Setup

        # CPO class
        self.setextra = setextra # set Eij, E_CAFFE, pfJ, etc. on every state update?
        super().__init__(mesh, boundaries, L, nu_multiplier=nu_multiplier, nu_realspace=nu_realspace, modelplane=modelplane, ds=ds, nvec=nvec)

        # specfab
        self.enhancementmodel = enhancementmodel # 'LTS', 'APEX', ...
        self.CAFFE_params = CAFFE_params # (Emin, Emax)
        self.n_EIE = n_EIE
        self.Lame_grain = self.sf.Cij_to_Lame_tranisotropic(Cij) 
        self.rho = rho

        # EnhancementFactor class
        self.symframe = symframe
        if   self.enhancementmodel == 'LTS':  homoparams = (Eij_grain, alpha, n_grain)
        elif self.enhancementmodel == 'APEX': homoparams = CAFFE_params
        else: homoparams = ()
        self.enhancementfactor = EnhancementFactor(mesh, L, enhancementmodel=enhancementmodel, homoparams=homoparams, modelplane=modelplane, symframe=symframe)

        ### Finish up
        
        self.initialize() # isotropic
        self.set_BCs([], [], []) # no BCs

        
    def initialize(self, *args, **kwargs):
        super().initialize(*args, **kwargs)
        self._setaux()
        
    def set_state(self, *args, **kwargs):
        super().set_state(*args, **kwargs)
        self._setaux()
        
    def evolve(self, u, *args, **kwargs):
        super().evolve(u, *args, **kwargs)
        self._setaux(u=u)

    def solvesteady(self, u, S, **kwargs):
        super().solvesteady(u, S, **kwargs)
        self._setaux(u=u)
        
    def _setaux(self, u=None):
        if self.setextra:
            self.mi, self.Eij, self.lami = self.get_Eij()
            self.xi, self.Exij, _        = self.get_Eij(ei=np.eye(3))
            # unpack above eigenframe fields for convenience
            self.m1, self.m2, self.m3 = self.mi # rheological symmetry directions
            self.E11, self.E22, self.E33, self.E23, self.E31, self.E12 = self.Eij  # eigenenhancements
            self.Exx, self.Eyy, self.Ezz, self.Eyz, self.Exz, self.Exy = self.Exij # Cartesian enhancements
            self.lam1, self.lam2, self.lam3 = self.lami # fabric eigenvalues \lambda_i (not a2 eigenvalues unless symframe=-1)
            # additional fields...
            self.pfJ = self.get_pfJ()
            if u is not None: 
                self.E_CAFFE = self.get_E_CAFFE(u)
                self.E_EIE   = self.get_E_EIE(u, self.Eij, self.mi)
            
    def get_pfJ(self, *args, **kwargs):
        return self.enhancementfactor.pfJ(self.s, *args, **kwargs)
            
    def get_E_CAFFE(self, u, *args, **kwargs):
        return self.enhancementfactor.E_CAFFE(self.s, u, Emin=self.CAFFE_params[0], Emax=self.CAFFE_params[1])

    def get_E_EIE(self, u, Eij, mi, *args, **kwargs):
        return self.enhancementfactor.E_EIE(u, Eij, mi, self.n_EIE)

    def get_Eij(self, ei=(), **kwargs):
        return self.enhancementfactor.Eij_tranisotropic(self.s, ei=ei, **kwargs)

    def get_elastic_velocities(self, x,y, theta,phi, alpha=1):
        nlm = self.get_state(x,y)
        vS1, vS2, vP = sf__.Vi_elastic_tranisotropic(nlm, alpha,self.Lame_grain,self.rho, theta,phi) # calculate elastic phase velocities using specfab
        return (vP, vS1, vS2)

    def Gamma0(self, u, T, A, Q):
        ### DDRX rate factor (generic)
        R = Constant(8.314) # gas constant (J/mol*K)
        D = sym(grad(u))
        epsE = sqrt(inner(D,D)/2)
        return project(epsE*Constant(A)*exp(-Constant(Q)/(R*T)), self.R)
        
    def Gamma0_Lilien23_EDC(self, u, T, A=4.3e7, Q=3.36e4):
        return self.Gamma0(u, T, A, Q) # EDC calibration by Lilien et al., 2023, p7

    def Gamma0_Lilien23_lab(self, u, T, A=1.91e7, Q=3.36e4):
        return self.Gamma0(u, T, A, Q) # Lab calibration by Lilien et al., 2023, p7
        
    def df2np(self, F, ele=('CG',2), withcoords=False):
        coords = copy.deepcopy(self.mesh.coordinates().reshape((-1, 2)).T)
        triang = tri.Triangulation(*coords, triangles=self.mesh.cells())    
        Q = FunctionSpace(self.mesh, *ele)
        F_np = project(F,Q).compute_vertex_values(self.mesh)
        return (triang, F_np, coords) if withcoords else (triang, F_np)
        
