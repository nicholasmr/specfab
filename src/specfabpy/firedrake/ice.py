#!/usr/bin/python3
# Nicholas Rathmann and Daniel Shapero, 2024

r"""
Firedrake all-in-one interface for ice fabric dynamics, viscous anisotropy, etc.

------
NOTICE
------
The 2D ice-flow model plane is assumed to be xz. 
This has the benifit of reducing the DOFs of the fabric evolution problem to involve only real numbers (real-valued state vector, s), considerably speeding up the solver.
If your problem is in fact xy, nothing changes in the way this class is used; the ODFs etc will simply reflect to assumption that flow is in xz rather than xy.
When visualizing ODFs/fabric quantities, you might therefore want to rotate the frame-of-reference from xz to xy.

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
"""

import numpy as np
from ..specfabpy import specfabpy as sf__
from .. import common as sfcom
from firedrake import *

class IceFabric:
    def __init__(
        self, mesh, boundaries, L, nu_multiplier=1, nu_realspace=1e-3, modelplane='xz', symframe=-1, ds=None, nvec=None
    ):
        ### Check args
        
        if modelplane != 'xz':
            raise ValueError('modelplane "%s" not supported, must be "xz"'%(modelplane))

        ### Setup
        
        self.mesh, self.boundaries = mesh, boundaries
        self.ds = ds   if ds   is not None else Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        self.n  = nvec if nvec is not None else FacetNormal(self.mesh)
        self.L = int(L) # spectral truncation
        self.symframe   = symframe
        self.modelplane = modelplane
        self.nu_realspace  = Constant(nu_realspace)
        self.nu_multiplier = Constant(nu_multiplier)

        ### Initialize fortran module
        
        self.sf = sf__
        self.lm, self.nlm_len = self.sf.init(self.L)
        self.rnlm_len = self.sf.get_rnlm_len()

        ### Fabric state and dyamics
        
        ele = ('CG',1)
        self.S = VectorFunctionSpace(self.mesh, *ele, dim=self.rnlm_len)
        self.R = FunctionSpace(self.mesh, *ele)
        self.G = TensorFunctionSpace(self.mesh, *ele, shape=(2,2))
        self.numdofs = self.R.dim()
        # Test, trial, solution container
        self.p = TrialFunction(self.S)
        self.q = TestFunction(self.S)
        self.ps = split(self.p)
        self.qs = split(self.q)
        self.s  = Function(self.S)
        self.s0 = Function(self.S)
        # Dynamical matrices rows, e.g. M_LROT[i,:] (this account for the real-real coefficient interactions, sufficient for xz problems)
        self.Mrr_LROT     = [Function(self.S) for ii in np.arange(self.rnlm_len)] # list of M matrix rows
        self.Mrr_DDRX_src = [Function(self.S) for ii in np.arange(self.rnlm_len)]
        self.Mrr_CDRX     = [Function(self.S) for ii in np.arange(self.rnlm_len)] 
        self.Mrr_REG      = [Function(self.S) for ii in np.arange(self.rnlm_len)] 

        ### Derived quantities: viscous anisotropy, a2 eigenvalues, J index, ...
        
        ele = ('DG',1)
        self.Sd = VectorFunctionSpace(self.mesh, *ele, dim=self.rnlm_len)
        self.Rd = FunctionSpace(self.mesh, *ele)
        self.Gd = TensorFunctionSpace(self.mesh, *ele, shape=(2,2))
        self.Vd = VectorFunctionSpace(self.mesh, *ele, dim=3) # for vectors
        self.numdofs0 = self.Rd.dim()
        self.mi   = [Function(self.Vd) for _ in range(3)] # (m1,m2,m3) fabric principal directions
        self.Eij  = [Function(self.Rd) for _ in range(6)] # Eij enhancement tensor
        self.lami = [Function(self.Rd) for _ in range(3)] # a2 eigenvalues (lami)
        self.E_CAFFE = Function(self.Rd)
        self.J = Function(self.Rd)

        ### Idealized states

        self.rnlm_iso  = [1/np.sqrt(4*np.pi)] + [0]*(self.rnlm_len-1) # Isotropic and normalized state
        self.nlm_dummy = np.zeros((self.nlm_len))
        self.M_zero    = np.zeros((self.rnlm_len,self.rnlm_len))

        ### Finish up
                
        self.initialize()
        
    def initialize(self, s0=None):
        s0_ = self.rnlm_iso if s0 is None else s0
        self.s.assign(project(Constant(s0_), self.S))
        self._nlm_nodal()

    def set_state(self, s):
        self.s.assign(s)
        self.s0.assign(s)

    def set_BCs(self, si, domids):
        self.bcs = [DirichletBC(self.S, si[ii], did) for ii, did in enumerate(domids)]

    def set_isotropic_BCs(self, domids):
        self.set_BCs([Constant(self.rnlm_iso)]*len(domids), domids)
        
    def get_nlm(self, *p, xz2xy=False):
        nlm = self.sf.rnlm_to_nlm(self.s(p)+0j, self.nlm_len)
        return self.sf.rotate_nlm_xz2xy(nlm) if xz2xy else nlm 
       
    def _nlm_nodal(self):
        sp = project(self.s, self.Sd)
        self.rnlm = np.array([sp.sub(ii).vector()[:] + 0j for ii in range(self.rnlm_len)]) # reduced form (rnlm) per node
        self.nlm  = np.array([self.sf.rnlm_to_nlm(self.rnlm[:,nn], self.nlm_len) for nn in range(self.numdofs0)]) # full form (nlm) per node (nlm[node,coef])

    def evolve(self, u, S, dt, iota=+1, Gamma0=None, Lambda0=None, steadystate=False):
        """
        Full-Stokes fabric solver
        """
        self.s0.assign(self.s)
        F = self._get_weakform(u, S, dt, iota, Gamma0, Lambda0, steadystate=steadystate)
        solve(lhs(F)==rhs(F), self.s, self.bcs, solver_parameters={'linear_solver':'gmres',}) # Non-symmetric system. Fastest tested are: gmres, bicgstab, tfqmr
        self._nlm_nodal()

    def _get_weakform(self, u, S, dt, iota, Gamma0, Lambda0, zeta=0, steadystate=False):

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
        D3f = np.array([self.mat3d(Df[nn]) for nn in np.arange(self.numdofs)])
        W3f = np.array([self.mat3d(Wf[nn]) for nn in np.arange(self.numdofs)])
        S3f = np.array([self.mat3d(Sf[nn]) for nn in np.arange(self.numdofs)])

        # Dynamical matrices M_* at each node as np arrays
        M_LROT     = np.array([self._M_reduced(self.sf.M_LROT,     self.nlm_dummy, D3f[nn], W3f[nn], iota, zeta) for nn in np.arange(self.numdofs)] ) if ENABLE_LROT else self.M_zero
        M_DDRX_src = np.array([self._M_reduced(self.sf.M_DDRX_src, self.nlm_dummy, S3f[nn])                      for nn in np.arange(self.numdofs)] ) if ENABLE_DDRX else self.M_zero
        M_REG      = np.array([self._M_reduced(self.sf.M_REG,      self.nlm_dummy, D3f[nn])                      for nn in np.arange(self.numdofs)] ) if ENABLE_REG  else self.M_zero

        # Populate entries of dynamical matrices row-wise
        kk = 0 # real-real interactions only in xz model frame
        for ii in np.arange(self.rnlm_len):
            if ENABLE_LROT: self.Mrr_LROT[ii].vector()[:]     = M_LROT[:,kk,ii,:]
            if ENABLE_DDRX: self.Mrr_DDRX_src[ii].vector()[:] = M_DDRX_src[:,kk,ii,:]
            if ENABLE_REG:  self.Mrr_REG[ii].vector()[:]      = M_REG[:,kk,ii,:]

        ### Construct weak form

        # Real space advection
        F = dot(dot(u, nabla_grad(self.p)), self.q)*dx

        # Time derivative
        dtinv = Constant(1/dt)
        if not steadystate:
            F += dtinv * dot( (self.p-self.s0), self.q)*dx

        # Real space stabilization (Laplacian diffusion)
        if ENABLE_REG:
            F += self.nu_realspace * inner(grad(self.p), grad(self.q))*dx

        # Lattice rotation
        if ENABLE_LROT:
            F += -sum([ dot(self.Mrr_LROT[ii], self.p)*self.qs[ii]*dx for ii in np.arange(self.rnlm_len)])

        # DDRX
        if ENABLE_DDRX:
            # @TODO NOT WORKING???
            F_src  = -Gamma0*sum([ dot(self.Mrr_DDRX_src[ii], self.p)*self.qs[ii]*dx for ii in np.arange(self.rnlm_len)]) 
            F_sink = -Gamma0*sum([ dot(self.Mrr_DDRX_src[0], self.s0)/self.s0[0]*self.ps[ii] * self.qs[ii]*dx for ii in np.arange(self.rnlm_len)]) # nonlinear sink term is linearized around previous solution (self.s0) following Rathmann and Lilien (2021)
            F += F_src - F_sink
#            F += F_sink

        # Orientation space stabilization (hyper diffusion)
        F += -self.nu_multiplier * sum([ dot(self.Mrr_REG[ii], self.p)*self.qs[ii]*dx for ii in np.arange(self.rnlm_len)])

        return F
        
    def _M_reduced(self, Mfun, *args, **kwargs):
        return self.sf.reduce_M(Mfun(*args, **kwargs), self.rnlm_len) # full form -> reduced form of M_*

    def mat3d(self, mat2d):
        return sfcom.mat3d(mat2d, self.modelplane, reshape=True)
        
    def eigenframe(self, *p, extrapolate=False):
        """
        Extract a2 eigenvectors (mi) and eigenvalues (lami) at specified point(s)
        
        Returns: (mi, lami)

        Note that you can back-construct a2 = lam1*np.outer(m1,m1) + lam2*np.outer(m2,m2) + lam3*np.outer(m3,m3) 
        """
        xf, yf = np.array(p[0], ndmin=1), np.array(p[1], ndmin=1)
        if len(xf.shape) > 1:
            raise ValueError('eigenframe(): please flatten input, only 1D (x,y) is supported')
        N = len(xf)
        mi, lami = np.zeros((N,3,3)), np.zeros((N,3))
        for ii in np.arange(N): 
            mi[ii], lami[ii] = sfcom.eigenframe(self.get_nlm(xf[ii],yf[ii]), symframe=self.symframe, modelplane=self.modelplane)
        return (mi[0], lami[0]) if N==1 else (mi, lami)
        
    def pfJ(self, *args, **kwargs):
        """
        Pole figure J (pfJ) index
        
        Returns: J
        """
        self.J.vector()[:] = sfcom.pfJ(self.nlm.T, *args, **kwargs)[:]
        return self.J
 
    def get_E_CAFFE(self, u, Emin=0.1, Emax=10):
        """
        CAFFE model (Placidi et al., 2010)
        
        Returns: E_CAFFE
        """
        Df = project(sym(grad(u)), self.Gd).vector()[:]
        self.E_CAFFE.vector()[:] = np.array([sf__.E_CAFFE(self.mat3d(Df[nn]), self.nlm[nn], Emin, Emax) for nn in np.arange(self.numdofs0)])
        return self.E_CAFFE
        
    def get_Eij(self, Eij_grain, alpha, n_grain, ei=()):   
        """
        Bulk enhancement factors w.r.t. ei=(e1,e2,e3) axes for *transversely isotropic* grains.
        If ei=() then CPO a2 eigenframe is used, ei=(m1,m2,m3) and returned Eij are the eigenenhancements.

        Returns: (mi, Eij, lami)
        """
        ### Check args
        
        if n_grain != 1: 
            raise ValueError('only n_grain = 1 (linear viscous) is supported')
        if not(0 <= alpha <= 1): 
            raise ValueError('alpha should be between 0 and 1')

        ### Calculate Eij etc. using specfabpy per node
        
        mi   = np.zeros((3, 3, self.numdofs0)) # (xyz component, i-th vector, node)
        Eij  = np.zeros((6, self.numdofs0))    # (Eij component, node)
        lami = np.zeros((3, self.numdofs0))    # (i-th a^(2) eigenvalue, node)        
        
        for nn in np.arange(self.numdofs0): 
            eigvecs, lami[:,nn] = sfcom.eigenframe(self.nlm[nn,:], symframe=self.symframe, modelplane=self.modelplane)
            mi[:,0,nn], mi[:,1,nn], mi[:,2,nn] = eigvecs.T if len(ei) == 0 else ei
            Eij[:,nn] = self.sf.Eij_tranisotropic(self.nlm[nn,:], mi[:,0,nn], mi[:,1,nn], mi[:,2,nn], Eij_grain, alpha, n_grain)
        
        """
        The enhancement-factor model depends on effective (homogenized) grain parameters, calibrated against deformation tests.
        For CPOs far from the calibration states, negative values *may* occur where Eij should tend to zero if truncation L is not large enough.
        """
        Eij[Eij<0] = 1e-2 # Set negative E_ij to a very small value (flow inhibiting)
        
        ### Populate functions
        
        for ii in range(3): # m1,m2,m3
            self.lami[ii].vector()[:] = lami[ii]
            for jj in range(3): # x,y,z
                self.mi[ii].sub(jj).vector()[:] = mi[jj,ii]
                
        for kk in range(6): 
            self.Eij[kk].vector()[:] = Eij[kk]
            
        return (self.mi, self.Eij, self.lami)
                
