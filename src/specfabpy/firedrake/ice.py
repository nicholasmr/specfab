#!/usr/bin/python3
# Nicholas Rathmann and Daniel Shapero, 2024

r"""
Firedrake all-in-one interface for ice fabric dynamics, viscous anisotropy, etc.
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
        self.symframe    = symframe
        self.modelplane  = modelplane
        # Real-space stabilization (multiplicative constant of real-space Laplacian)
        self.nu_realspace = Constant(nu_realspace)
        # Multiplier of orientation-space regularization magnitude
        self.nu_multiplier = Constant(nu_multiplier)

        ### Initialize specfab fortran module
        self.sf = sf__
        self.lm, self.nlm_len_full = self.sf.init(self.L)
        self.USE_REDUCED = True 
        self.nlm_len = self.sf.get_rnlm_len() if self.USE_REDUCED else self.nlm_len_full

        ### Fabric dyamics
        eletype, eleorder = 'CG', 1
        self.Rele = FiniteElement(eletype, self.mesh.ufl_cell(), eleorder)
        self.Vele = MixedElement([self.Rele]*self.nlm_len)
        self.R = FunctionSpace(self.mesh, self.Rele)
        self.V = VectorFunctionSpace(self.mesh, "CG", 1, dim=self.nlm_len)
        self.G = TensorFunctionSpace(self.mesh, eletype, eleorder, shape=(2,2))
        # Get vector() size w/ MPI (else self.R.dim())
        self.numdofs = Function(self.R).vector().local_size()

        ### Viscous anisotropy
        self.R0 = FunctionSpace(self.mesh, "DG", 0)
        self.G0 = TensorFunctionSpace(self.mesh, "DG", 0, shape=(2,2))
        self.numdofs0 = Function(self.R0).vector().local_size()

        ### Function spaces
        self.Wele = self.Vele
        self.W = self.V
        self.pr = TrialFunction(self.W) # unknowns
        self.qr = TestFunction(self.W)  # weight functions
        self.pr_sub = split(self.pr)    # easy access of subelements
        self.qr_sub = split(self.qr)    # easy access of subelements
        # Solution containers
        self.w      = Function(self.W)
        self.w_prev = Function(self.W)
        # Dynamical matrices (real-real interactions only for xz problem)
        self.Mrr_LROT     = [Function(self.V) for ii in np.arange(self.nlm_len)]
        self.Mrr_DDRX_src = [Function(self.V) for ii in np.arange(self.nlm_len)]
        self.Mrr_CDRX     = [Function(self.V) for ii in np.arange(self.nlm_len)] 
        self.Mrr_REG      = [Function(self.V) for ii in np.arange(self.nlm_len)] 

        ### Idealized states
        self.nlm_iso   = [1/np.sqrt(4*np.pi)] + [0]*(self.nlm_len-1) # Isotropic and normalized state
        self.nlm_zero  = [0]*(self.nlm_len)
        self.nlm_dummy = np.zeros((self.nlm_len_full))
        
    def initialize(self, wr=None, wi=None):
        """
        Initialize uniform CPO field
        """
        if wr is None: wr, wi = self.nlm_iso, self.nlm_zero
        self.w.assign(project(Constant(wr), self.V)) # real part

    def set_state(self, w):
        self.w.assign(w)
        self.w_prev.assign(w)

    def set_BCs(self, wr, wi, domids):
        """
        Set boundary conditions
        """
        self.bcs = [DirichletBC(self.W, wr[ii], did) for ii, did in enumerate(domids)] # real part

    def set_isotropic_BCs(self, domids):
        """
        Easy setting of isotropic boundary conditions
        """
        wr = [Constant(self.nlm_iso)]  * len(domids)
        wi = [Constant(self.nlm_zero)] * len(domids)
        self.set_BCs(wr, wi, domids)

    def weakform(self, u, S, dt, iota, Gamma0, Lambda0, zeta=0, steadystate=False):
        """
        Build weak form from dynamical matrices
        """

        ENABLE_LROT = iota is not None
        ENABLE_DDRX = (Gamma0 is not None) and (S is not None)
        ENABLE_CDRX = Lambda0 is not None

        # Flattened strain-rate and spin tensors for accessing them per node
        Df = project( sym(grad(u)), self.G).vector()[:] # strain rate
        Wf = project(skew(grad(u)), self.G).vector()[:] # spin
        Sf = project(S, self.G).vector()[:] # deviatoric stress

        # Dynamical matrices at each DOF (note indexing is different from FEniCS)
        if ENABLE_LROT:
            M_LROT_nodal = np.array([self.sf.reduce_M(self.sf.M_LROT(self.nlm_dummy, self.mat3d(Df[nn]), self.mat3d(Wf[nn]), iota, zeta), self.nlm_len) for nn in np.arange(self.numdofs)] )
        if ENABLE_DDRX: 
            M_DDRX_src_nodal = np.array([self.sf.reduce_M(self.sf.M_DDRX_src(self.nlm_dummy, self.mat3d(Sf[nn])), self.nlm_len) for nn in np.arange(self.numdofs)] )
        M_REG_nodal = np.array([self.sf.reduce_M(self.sf.M_REG(self.nlm_dummy, self.mat3d(Df[nn])), self.nlm_len) for nn in np.arange(self.numdofs)] )

        # Populate entries of dynamical matrices
        kk = 0 # real-real interactions only
        for ii in np.arange(self.nlm_len):
            if ENABLE_LROT: self.Mrr_LROT[ii].vector()[:] = M_LROT_nodal[:,kk,ii,:]
            if ENABLE_DDRX: self.Mrr_DDRX_src[ii].vector()[:] = M_DDRX_src_nodal[:,kk,ii,:]
            self.Mrr_REG[ii].vector()[:] = M_REG_nodal[:,kk,ii,:]

        ### Construct weak form

        # Real space advection, div(s*u)
        F = dot(dot(u, nabla_grad(self.pr)), self.qr)*dx # real part

        # Time derivative
        dtinv = Constant(1/dt)
        if not steadystate:
            F += dtinv * dot( (self.pr-self.w_prev), self.qr)*dx # real part

        # Real space stabilization (Laplacian diffusion)
        F += self.nu_realspace * inner(grad(self.pr), grad(self.qr))*dx # real part

        # Lattice rotation
        if ENABLE_LROT:
            F += -sum([ dot(self.Mrr_LROT[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # real part

        # DDRX
        if ENABLE_DDRX:
            F_src  = sum([ -Gamma0*dot(self.Mrr_DDRX_src[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)])
            # nonlinear sink term is linearized around previous solution (self.w_prev) following Rathmann and Lilien (2021)
            F_sink = sum([ -Gamma0*dot(Mrr_DDRX_src[0],self.w_prev)/self.w_prev[0]*self.pr_sub[ii] * self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)])
            #F_sink = sum([ -Gamma0*dot(self.Mrr_DDRX_src[0],self.pr)*self.pr_sub[ii] * self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # nonlinear sink term is linearized around previous solution (self.w_prev) following Rathmann and Lilien (2021)
            F += F_src - F_sink

        # Orientation space stabilization (hyper diffusion)
        F += -self.nu_multiplier * sum([ dot(self.Mrr_REG[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # real part

        return F

    def evolve(self, u, S, dt, iota=+1, Gamma0=None, Lambda0=None, steadystate=False):
        """
        Evolve CPO using Laplacian stabilized, Euler time integration
        """
        if self.w is None:
            raise ValueError('CPO state "w" not set. Did you forget to initialize the CPO field?')
        self.w_prev.assign(self.w) # current state (w) must be set
        F = self.weakform(u, S, dt, iota, Gamma0, Lambda0, steadystate=steadystate)
        solve(lhs(F)==rhs(F), self.w, self.bcs, solver_parameters={'linear_solver':'gmres', }) # fastest tested are: gmres, bicgstab, tfqmr (For non-symmetric problems, a Krylov solver for non-symmetric systems, such as GMRES, is a better choice)
        
    def get_nlm(self, x, y):    
        """
        Extract full-form CPO state vector at point
        """        
        nlm = self.w((x,y)) + 0j
        return self.sf.rnlm_to_nlm(nlm, self.nlm_len_full) if self.USE_REDUCED else nlm
       
    def eigenframe(self, x, y, extrapolate=False):
        """
        Extract a2 eigenvectors and eigenvalues at point
        @TODO rewrite so that it returns "FunctionSpace variables"
        """
#        self.w.set_allow_extrapolation(extrapolate) # @TODO not supported in firedrake?
        xf, yf = np.array(x, ndmin=1), np.array(y, ndmin=1)
        if len(xf.shape) > 1: 
            raise ValueError('eigenframe(): only 1D (x,y) are supported (i.e. flattened)')
        N = len(xf)
        eigvecs, eigvals = np.zeros((N,3,3)), np.zeros((N,3))
        for ii in np.arange(N): 
            eigvecs[ii,:,:], eigvals[ii,:] = sfcom.eigenframe(self.get_nlm(xf[ii],yf[ii]), symframe=self.symframe, modelplane=self.modelplane)
        return (eigvecs[0,:,:], eigvals[0,:]) if N==1 else (eigvecs, eigvals)
    
    def mat3d(self, mat2d):
        return sfcom.mat3d(mat2d, self.modelplane, reshape=True) # from common.py
    
    def nlm_nodal(self):
        """
        state vector nlm per node as np array
        """
        rnlm = np.array([ project(self.w[ii],self.R0).vector()[:] + 0j for ii in range(self.nlm_len) ]) # reduced form per node
        nlm = np.array([ self.sf.rnlm_to_nlm(rnlm[:,nn], self.nlm_len_full) for nn in range(self.numdofs0) ]) # *full* form
        return nlm # nlm[node,coef]
 
    def E_CAFFE(self, u, Emin=0.1, Emax=10):
        """
        CAFFE model (Placidi et al., 2010)
        """
        Df = project( sym(grad(u)), self.G0).vector()[:] # flattened strain-rate tensor
        nlm = self.nlm_nodal() # (nlm component, node)        
        E_df = Function(self.R0)
        E_df.vector()[:] = np.array([sf__.E_CAFFE(self.mat3d(Df[nn]), nlm[nn], Emin, Emax) for nn in np.arange(self.numdofs0) ])
        return E_df
        
