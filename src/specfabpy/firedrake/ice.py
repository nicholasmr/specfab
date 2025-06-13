#!/usr/bin/python3
# Nicholas Rathmann and Daniel Shapero, 2024

"""
Firedrake all-in-one interface for ice fabric dynamics, viscous anisotropy, etc.

-----------
NOTICE
-----------
The 2D ice-flow model plane is assumed to be xz. 
This has the benifit of reducing the DOFs of the fabric evolution problem to involve only real numbers (real-valued state vector, s), considerably speeding up the solver.
If your problem is in fact xy, nothing changes in the way this class is used; the (M)ODFs etc will simply reflect to assumption that flow is in xz rather than xy.
When visualizing (M)ODFs/fabric quantities, you might therefore want to rotate the frame-of-reference from xz to xy.

-----------
DEFINITIONS
-----------
s               : fabric state vector [n_0^0, n_2^0, n_2^1, n_2^2, ...], where n_l^m (nlm in code) are the harmonic expansion coefficients
a2              : 2nd order orientation (structure) tensor 
lami            : eigenvlues of a2
mi              : fabric principal directions {m1,m2,m3} (eigenvectors of a2)
Eij             : diretional enhancement factors w.r.t. some frame {e1,e2,e3}, e.g. {m1,m2,m3}
nu_realspace    : real-space stabilization (multiplicative constant of real-space Laplacian)
nu_multiplier   : multiplier of orientation-space regularization magnitude
iota            : plastic spin free parameter; iota=1 => deck-of-cards behaviour
Gamma0          : DDRX rate factor so that Gamma=Gamma0*(D-<D>), where D is the deformability
Lambda0         : CDRX rate factor

-----------
TODO
-----------
- SSA fabric source/sinks and topo advection terms need to be included 
- ...
"""

import numpy as np
import matplotlib.tri as tri
#import code # code.interact(local=locals())

from ..specfabpy import specfabpy as sf__
from .. import constants as sfconst
from .. import common as sfcom

from firedrake import *

class IceFabric:
    def __init__(
        self, mesh, boundaries, L, *args, 
        nu_multiplier=1, nu_realspace=1e-3, modelplane='xz', symframe=-1, 
        Eij_grain=(1,1), alpha=0, n_grain=1, CAFFE_params=(0.1, 10, 1), n_EIE=3, 
        ds=None, nvec=None, setextra=True, 
        Cij=sfconst.ice['elastic']['Bennett1968'], rho=sfconst.ice['density'], **kwargs
    ):
        ### Check args
        
        if modelplane != 'xz':
            raise ValueError('modelplane "%s" not supported, must be "xz"'%(modelplane))
            
        if n_grain != 1:
            raise ValueError('only n_grain = 1 (linear viscous) is supported')
            
        if not(0 <= alpha <= 1):
            raise ValueError('alpha should be between 0 and 1')

        ### Setup
        
        # Firedrake
        self.mesh, self.boundaries = mesh, boundaries
        self.ds = ds   if ds   is not None else Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        self.n  = nvec if nvec is not None else FacetNormal(self.mesh)
        self.setextra = setextra # set Eij, E_CAFFE, pfJ, etc. on every state update?

        # specfab
        self.L = int(L) # spectral truncation
        self.sf = sf__
        self.lm, self.nlm_len = self.sf.init(self.L)
        self.rnlm_len = self.sf.get_rnlm_len()
        self.symframe   = symframe
        self.modelplane = modelplane
        self.nu_realspace  = nu_realspace
        self.nu_multiplier = nu_multiplier
        self.grain_params = (Eij_grain, alpha, n_grain)
        self.CAFFE_params = CAFFE_params # (Emin, Emax)
        self.n_EIE        = n_EIE
        self.Lame_grain   = self.sf.Cij_to_Lame_tranisotropic(Cij) 
        self.rho          = rho

        ### Fabric state and dyamics
        
        ele = ('CG',1)
        self.S = VectorFunctionSpace(self.mesh, *ele, dim=self.rnlm_len)
        self.R = FunctionSpace(self.mesh, *ele)
        self.G = TensorFunctionSpace(self.mesh, *ele, shape=(2,2))
        self.numdofs = self.R.dim()
        
        # Test, trial, solution container
        self.p = TrialFunction(self.S)
        self.w = TestFunction(self.S)
        self.ps = split(self.p)
        self.ws = split(self.w)
        self.s  = Function(self.S)
        self.s0 = Function(self.S)
        
        # Dynamical matrices rows, e.g. M_LROT[i,:] (this account for the real-real coefficient interactions, sufficient for xz problems)
        self.Mrr_LROT     = [Function(self.S) for ii in range(self.rnlm_len)] # list of M matrix rows
        self.Mrr_DDRX_src = [Function(self.S) for ii in range(self.rnlm_len)]
        self.Mrr_CDRX     = [Function(self.S) for ii in range(self.rnlm_len)] 
        self.Mrr_REG      = [Function(self.S) for ii in range(self.rnlm_len)] 

        ### Derived quantities: viscous anisotropy, a2 eigenvalues, J index, ...
        
        ele = ('DG',0) # safe for forward modeling of coupled fabric--flow problems
        self.Sd = VectorFunctionSpace(self.mesh, *ele, dim=self.rnlm_len)
        self.Rd = FunctionSpace(self.mesh, *ele)
        self.Gd = TensorFunctionSpace(self.mesh, *ele, shape=(2,2))
        self.Vd = VectorFunctionSpace(self.mesh, *ele, dim=3) # for vectors
        self.numdofs0 = self.Rd.dim()

        ### Idealized states

        self.rnlm_iso  = [1/np.sqrt(4*np.pi)] + [0]*(self.rnlm_len-1) # Isotropic and normalized state
        self.nlm_dummy = np.zeros((self.nlm_len))
        self.M_zero    = np.zeros((self.numdofs, self.rnlm_len, self.rnlm_len))

        ### Aux/shortcuts
        
        self.dofs  = np.arange(self.numdofs)
        self.dofs0 = np.arange(self.numdofs0)
        self.srng  = np.arange(self.nlm_len)
        self.srrng = np.arange(self.rnlm_len)

        ### Finish up
                
        self.initialize()


    def initialize(self, s0=None):
        s0 = self.rnlm_iso if s0 is None else s0
        self.set_state(project(Constant(s0), self.S))

    def set_state(self, s):
        self.s.assign(s)
        self.s0.assign(s)
        self._setaux()
        
    def set_BCs(self, si, domids):
        self.bcs = [DirichletBC(self.S, si[ii], did) for ii, did in enumerate(domids)]

    def set_isotropic_BCs(self, domids):
        self.set_BCs([Constant(self.rnlm_iso)]*len(domids), domids)
        
    def evolve(self, u, S, dt, iota=+1, Gamma0=None, Lambda0=None):
        """
        The fabric solver, called to step fabric field forward in time by amount dt
        """
        self.s0.assign(self.s)
        F = self._weakform(self.p, u,S, dt, iota,Gamma0,Lambda0)
        solve(lhs(F)==rhs(F), self.s, self.bcs, solver_parameters={'linear_solver':'gmres',}) # non-symmetric system (fastest tested are: gmres, bicgstab, tfqmr)
        self._setaux(u=u)
        
    def solvesteady(self, u, S, iota=+1, Gamma0=None, Lambda0=None, LROT_guess=False, **kwargs):

        dt = None # select steady weak form

        # LROT solution (+CDRX+REG) *or* used as init guess for nonlinear problem
        if Gamma0 is None or LROT_guess:
            print('*** Linear solve for LROT-only problem %s'%('(used as initial guess)' if (LROT_guess and (Gamma0 is not None)) else ''))
            self.evolve(u, S, dt, iota=iota, Gamma0=None, Lambda0=Lambda0, **kwargs) # solution stored in self.s 

        # LROT + DDRX solution (+CDRX+REG)
        if Gamma0 is not None: # nonlinear problem?
            print('*** Nonlinear solve for DDRX-activated problem')
            if not isinstance(Gamma0, list): Gamma0 = [Gamma0,]
            for ii, gam0 in enumerate(Gamma0):
                print('...approaching solution using DDRX rate factor %i of %i'%(ii, len(Gamma0)))
                F = self._weakform(self.s, u, S, dt, iota,gam0,Lambda0, **kwargs)
                solve(F==0, self.s, self.bcs, solver_parameters={'newton_solver':{'linear_solver':'gmres', 'preconditioner':'none'}})

        self._setaux(u=u)

    def _weakform(self, s, u,S, dt, iota=None, Gamma0=None, Lambda0=None, zeta=0):

        # Crystal processes to include
        ENABLE_LROT = iota is not None
        ENABLE_DDRX = (Gamma0 is not None) and (S is not None)
        ENABLE_CDRX = Lambda0 is not None
        ENABLE_REG  = True
        
        # Flattened strain-rate and spin tensors for accessing them per node
        Df = project( sym(grad(u)), self.G).vector()[:] # strain rate
        Wf = project(skew(grad(u)), self.G).vector()[:] # spin
        Sf = project(S, self.G).vector()[:] # deviatoric stress

        # Same but in 3D for fabric problem
        D3 = self.mat3d(Df) # [node,3,3]
        W3 = self.mat3d(Wf)
        S3 = self.mat3d(Sf)

        # Dynamical matrices M_* at each node as np arrays (indexing is M[node,row,column])
        M_LROT     = np.array([self._M_reduced(self.sf.M_LROT,     self.nlm_dummy, D3[nn], W3[nn], iota, zeta) for nn in self.dofs] ) if ENABLE_LROT else self.M_zero
        M_DDRX_src = np.array([self._M_reduced(self.sf.M_DDRX_src, self.nlm_dummy, S3[nn])                     for nn in self.dofs] ) if ENABLE_DDRX else self.M_zero
        M_REG      = np.array([self._M_reduced(self.sf.M_REG,      self.nlm_dummy, D3[nn])                     for nn in self.dofs] ) if ENABLE_REG  else self.M_zero

        # Populate entries of dynamical matrices *row-wise* (row index is ii)
        for ii in self.srrng:
            if ENABLE_LROT: self.Mrr_LROT[ii].vector()[:]     = M_LROT[:,ii,:] # all nodes, row=ii, all columns
            if ENABLE_DDRX: self.Mrr_DDRX_src[ii].vector()[:] = M_DDRX_src[:,ii,:]
            if ENABLE_REG:  self.Mrr_REG[ii].vector()[:]      = M_REG[:,ii,:]

        ### Construct weak form

        # dummy zero term to make rhs(F) work when solving steady-state problem 
        # this can probably be removed once the SSA source/sink terms are added
        s_null = Function(self.S)
        s_null.vector()[:] = 0.0
        F = dot(s_null, self.w)*dx

        # Time derivative
        if dt is not None:
            F += Constant(1/dt) * dot( (s-self.s0), self.w)*dx

        # Real space advection
        F += dot(dot(u, nabla_grad(s)), self.w)*dx

        # Real space stabilization (Laplacian diffusion)
        F += Constant(self.nu_realspace) * inner(grad(s), grad(self.w))*dx

        if ENABLE_LROT:
            F += -sum([ dot(self.Mrr_LROT[ii], s)*self.ws[ii]*dx for ii in self.srrng])

        if ENABLE_DDRX:
            F_src  = -sum([Gamma0*dot(self.Mrr_DDRX_src[ii], s)*self.ws[ii]*dx for ii in self.srrng]) 
            
            s_sub = split(s) if dt is None else self.ps
            if dt is None: 
                # Steady state solver? Use nonlinear iteration for DDRX-activated problem
                D = dot(self.Mrr_DDRX_src[0], s)/s_sub[0] # <D> (Rathmann et al., 2025)                    
            else: 
                # Time-dependent probem? Linearized problem by linearizing <D>
                D = dot(self.Mrr_DDRX_src[0], self.s0)/self.s0.sub(0) # <D> estimated using previous solution (Rathmann and Lilien, 2021)
            F_snk = -sum([Gamma0*D*s_sub[ii]*self.ws[ii]*dx for ii in self.srrng])
            
            F += F_src - F_snk

        # Orientation space stabilization (hyper diffusion)
        if ENABLE_REG:
            F += -Constant(self.nu_multiplier)*sum([dot(self.Mrr_REG[ii], s)*self.ws[ii]*dx for ii in self.srrng])

        return F
        
    def _M_reduced(self, Mfun, *args, **kwargs):
        return self.sf.reduce_M(Mfun(*args, **kwargs), self.rnlm_len)[0] # full form -> reduced form of M_*

    def mat3d(self, D2):
        return sfcom.mat3d_arr(D2, self.modelplane, reshape=True) # can take array of D2
        
    def eigenframe(self, *p, **kwargs):
        """
        Extract a2 eigenvectors (mi) and eigenvalues (lami) at specified list of points p=([x1,x2,...], [y1,y2,...])
        Note that you can back-construct a2 = lam1*np.outer(m1,m1) + lam2*np.outer(m2,m2) + lam3*np.outer(m3,m3) 
                
        Returns: (mi, lami)
        """
        xf, yf = np.array(p[0], ndmin=1), np.array(p[1], ndmin=1) # (point num, x/y value)
        nlm = np.array([self.get_nlm(xf[pp],yf[pp], **kwargs) for pp in range(len(xf))]) # (point num, nlm coefs)
        mi, lami = sfcom.eigenframe(nlm, symframe=self.symframe, modelplane=self.modelplane)
        return (mi, lami)
        
    def get_nlm(self, *p, xz2xy=False):
        nlm = self.sf.rnlm_to_nlm(self.s(p)+0j, self.nlm_len)
        return self.sf.rotate_nlm_xz2xy(nlm) if xz2xy else nlm 
       
    def _setaux(self, u=None):
        """
        Set auxiliary fields
        """
        sp = project(self.s, self.Sd)
        self.rnlm = np.array([sp.sub(ii).vector()[:] + 0j for ii in self.srrng]) # reduced form (rnlm) per node
        self.nlm  = np.array([self.sf.rnlm_to_nlm(self.rnlm[:,nn], self.nlm_len) for nn in self.dofs0]) # full form (nlm) per node (nlm[node,coef])
        
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
                self.E_EIE   = self.get_E_EIE(u)
        
    def get_pfJ(self, *args, **kwargs):
        """
        Pole figure J (pfJ) index
        """
        pfJ = Function(self.Rd)
        pfJ.vector()[:] = sfcom.pfJ(self.nlm, *args, **kwargs)[:]
        return pfJ
 
    def get_E_CAFFE(self, u):
        """
        CAFFE model (Placidi et al., 2010)
        """
        Df = project(sym(grad(u)), self.Gd).vector()[:]
        E_CAFFE = Function(self.Rd)
        E_CAFFE.vector()[:] = self.sf.E_CAFFE_arr(self.nlm, self.mat3d(Df), *self.CAFFE_params)
        return E_CAFFE
        
    def get_E_EIE(self, u):
        """
        EIE model (Rathmann et al., in prep)
        *** Not yet implemented ***
        """
#        Df = project(sym(grad(u)), self.Gd).vector()[:]
        E_EIE = Function(self.Rd)
#        E_EIE.vector()[:] = self.sf.E_EIE_arr(...)
        return E_EIE
        
    def get_Eij(self, ei=()):   
        """
        Bulk enhancement factors wrt ei=(e1,e2,e3) axes for *transversely isotropic* grains.
        If ei=() then CPO eigenframe is used.
        """

        ### Calculate Eij etc. using specfabpy per node
        
        mi, lami = sfcom.eigenframe(self.nlm, symframe=self.symframe, modelplane=self.modelplane) # [node,i,xyz], [node,i]
        if   len(ei) == 0:           ei = (mi[:,0], mi[:,1], mi[:,2])
        elif len(ei[0].shape) == 1:  ei = sfcom.ei_tile(ei, self.numdofs0) # ei[i][node,xyz] = (e1[node,xyz], e2, e3)
        Eij = self.sf.Eij_tranisotropic_arr(self.nlm, *ei, *self.grain_params) # [node, Voigt index 1--6]

        # The enhancement factor model depends on effective (homogenized) grain parameters, calibrated against deformation tests.
        # For CPOs far from the calibration states, negative values *may* occur where Eij should tend to zero if truncation L is not large enough.
        Eij[Eij<0] = 1e-2 # Set negative E_ij to a very small value (flow inhibiting)
        
        ### Populate functions
        
        ei_fs   = [Function(self.Vd) for _ in range(3)] # a2 eigenvectors if ei=()
        lami_fs = [Function(self.Rd) for _ in range(3)] # a2 eigenvalues (lami)
        Eij_fs  = [Function(self.Rd) for _ in range(6)] # Eij enhancement tensor
        
        for ii in range(3): # mi=m1,m2,m3
            lami_fs[ii].vector()[:] = lami[:,ii]
            for jj in range(3): # j=x,y,z
                ei_fs[ii].sub(jj).vector()[:] = ei[ii][:,jj]
                
        for kk in range(6): 
            Eij_fs[kk].vector()[:] = Eij[:,kk]
            
        return (ei_fs, Eij_fs, lami_fs)
                
    def Gamma0_Lilien23_EDC(self, u, T, A=4.3e7, Q=3.36e4):
        # DDRX rate factor from Dome C ice-core calibration experiment (Lilien et al., 2023, p. 7)
        R = Constant(8.314) # gas constant (J/mol*K)
        D = sym(grad(u))
        epsE = sqrt(inner(D,D)/2)
        return project(epsE*Constant(A)*exp(-Constant(Q)/(R*T)), self.R)

    def Gamma0_Lilien23_lab(self, u, T, A=1.91e7, Q=3.36e4):
        # DDRX rate factor from lab calibration experiment (Lilien et al., 2023, p. 7)
        return self.Gamma0_Lilien23_EDC(u, T, A=A, Q=Q)

