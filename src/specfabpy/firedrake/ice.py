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
        self.symframe   = symframe
        self.modelplane = modelplane
        self.nu_realspace  = Constant(nu_realspace)  # Real-space stabilization (multiplicative constant of real-space Laplacian)
        self.nu_multiplier = Constant(nu_multiplier) # Multiplier of orientation-space regularization magnitude

        ### Initialize specfab fortran module
        self.sf = sf__
        self.lm, self.nlm_len_full = self.sf.init(self.L)
        self.USE_REDUCED = True
        self.nlm_len = self.sf.get_rnlm_len() if self.USE_REDUCED else self.nlm_len_full

        ### Fabric dyamics
        ele = ('CG', 1)
        self.R = FunctionSpace(self.mesh, *ele)
        self.W = VectorFunctionSpace(self.mesh, *ele, dim=self.nlm_len)
        self.G = TensorFunctionSpace(self.mesh, *ele, shape=(2,2))
        self.numdofs = Function(self.R).vector().local_size() # Get vector() size w/ MPI (else self.R.dim() is sufficient)

        ### Viscous anisotropy
#        ele = ('DG', 0)
        ele = ('CG', 1)
        self.R0  = FunctionSpace(self.mesh, *ele)
        self.R30 = VectorFunctionSpace(self.mesh, *ele, dim=3) # for vectors
        self.G0  = TensorFunctionSpace(self.mesh, *ele, shape=(2,2))
        self.numdofs0 = Function(self.R0).vector().local_size()

        ### Function spaces
        self.pr = TrialFunction(self.W) # unknowns
        self.qr = TestFunction(self.W)  # weight functions
        self.pr_sub = split(self.pr)    # easy access of subelements
        self.qr_sub = split(self.qr)    # easy access of subelements
        # Solution containers
        self.w      = Function(self.W)
        self.w_prev = Function(self.W)
        # Dynamical matrices (real-real interactions only for xz problem)
        self.Mrr_LROT     = [Function(self.W) for ii in np.arange(self.nlm_len)] # list of M matrix rows
        self.Mrr_DDRX_src = [Function(self.W) for ii in np.arange(self.nlm_len)]
        self.Mrr_CDRX     = [Function(self.W) for ii in np.arange(self.nlm_len)] 
        self.Mrr_REG      = [Function(self.W) for ii in np.arange(self.nlm_len)] 

        ### Idealized states
        self.nlm_iso   = [1/np.sqrt(4*np.pi)] + [0]*(self.nlm_len-1) # Isotropic and normalized state
        self.nlm_zero  = [0]*(self.nlm_len)
        self.nlm_dummy = np.zeros((self.nlm_len_full))
        
        self.initialize()
        
    def initialize(self, w0=None):
        """
        Initialize uniform CPO field
        """
        w0_ = self.nlm_iso if w0 is None else w0
        self.w.assign(project(Constant(w0_), self.W)) # real part

    def set_state(self, w):
        """
        Set CPO field state
        """
        self.w.assign(w)
        self.w_prev.assign(w)

    def set_BCs(self, w, domids):
        """
        Set boundary conditions
        """
        self.bcs = [DirichletBC(self.W, w[ii], did) for ii, did in enumerate(domids)] # real part

    def set_isotropic_BCs(self, domids):
        """
        Easy setting of isotropic boundary conditions
        """
        self.set_BCs([Constant(self.nlm_iso)]*len(domids), domids)

    def evolve(self, u, S, dt, iota=+1, Gamma0=None, Lambda0=None, steadystate=False):
        """
        Evolve CPO using Laplacian stabilized, Euler time integration
        """
        self.w_prev.assign(self.w) # current state (w) must be set
        F = self.weakform(u, S, dt, iota, Gamma0, Lambda0, steadystate=steadystate)
        solve(lhs(F)==rhs(F), self.w, self.bcs, solver_parameters={'linear_solver':'gmres',}) # fastest tested are: gmres, bicgstab, tfqmr --- note this is a non-symmetric system!

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
        kk = 0 # real-real interactions only in xz model frame
        for ii in np.arange(self.nlm_len):
            if ENABLE_LROT: self.Mrr_LROT[ii].vector()[:] = M_LROT_nodal[:,kk,ii,:]
            if ENABLE_DDRX: self.Mrr_DDRX_src[ii].vector()[:] = M_DDRX_src_nodal[:,kk,ii,:]
            self.Mrr_REG[ii].vector()[:] = M_REG_nodal[:,kk,ii,:]

        ### Construct weak form

        # Real space advection
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
            F_src  = sum([ -Gamma0*dot(Mrr_DDRX_src[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) 
            F_sink = sum([ -Gamma0*dot(Mrr_DDRX_src[0],self.w_prev)/self.w_prev[0]*self.pr_sub[ii] * self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # nonlinear sink term is linearized around previous solution (self.w_prev) following Rathmann and Lilien (2021)
            F += F_src - F_sink

        # Orientation space stabilization (hyper diffusion)
        F += -self.nu_multiplier * sum([ dot(self.Mrr_REG[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # real part

        return F

    def get_nlm(self, x, y):    
        """
        Extract full-form CPO state vector at point
        """        
        nlm = self.w((x,y)) + 0j
        return self.sf.rnlm_to_nlm(nlm, self.nlm_len_full) if self.USE_REDUCED else nlm
       
    def mat3d(self, mat2d):
        return sfcom.mat3d(mat2d, self.modelplane, reshape=True) # from common.py
    
    def nlm_nodal(self):
        """
        state vector "s" (i.e. coefficient n_l^m) per node as np array
        """
        rnlm = np.array([ project(self.w[ii],self.R0).vector()[:] + 0j for ii in range(self.nlm_len) ]) # reduced form per node
        nlm = np.array([ self.sf.rnlm_to_nlm(rnlm[:,nn], self.nlm_len_full) for nn in range(self.numdofs0) ]) # *full* form
        return nlm # nlm[node,coef]
        
    def eigenframe(self, x, y, extrapolate=False):
        """
        Extract a2 eigenvectors and eigenvalues (lami) at specified point(s)
        """
        xf, yf = np.array(x, ndmin=1), np.array(y, ndmin=1)
        if len(xf.shape) > 1: 
            raise ValueError('eigenframe(): please flatten input, only 1D (x,y) is supported')
        N = len(xf)
        mi, lami = np.zeros((N,3,3)), np.zeros((N,3))
        for ii in np.arange(N): 
            mi[ii,:,:], lami[ii,:] = sfcom.eigenframe(self.get_nlm(xf[ii],yf[ii]), symframe=self.symframe, modelplane=self.modelplane)
        return (mi[0,:,:], lami[0,:]) if N==1 else (mi, lami)
        
    def pfJ(self, **kwargs):
        """
        Pole figure J (pfJ) index
        """
        nlm = self.nlm_nodal() # (nlm component, node)
        J = pfJ(nlm, **kwargs) # sfcom.pfJ()
        J_df = Function(self.R0)
        J_df.vector()[:] = J[:] # nonzero imaginary parts are numerical uncertainty, should be real-valued
        return J_df
 
    def E_CAFFE(self, u, Emin=0.1, Emax=10):
        """
        CAFFE model (Placidi et al., 2010)
        """
        Df = project( sym(grad(u)), self.G0).vector()[:] # flattened strain-rate tensor
        nlm = self.nlm_nodal() # (nlm component, node)        
        E_df = Function(self.R0)
        E_df.vector()[:] = np.array([sf__.E_CAFFE(self.mat3d(Df[nn]), nlm[nn], Emin, Emax) for nn in np.arange(self.numdofs0) ])
        return E_df
        
    def Eij(self, Eij_grain, alpha, n_grain, ei_arg=()):   
        """
        Bulk enhancement factors w.r.t. ei=(e1,e2,e3) axes for *transversely isotropic* grains.
        If mi=() then CPO a2 eigenframe is used, ei=(m1,m2,m3) and returned Eij are the eigenenhancements.
        """
        ### Check args
        if n_grain != 1: 
            raise ValueError('only n_grain = 1 (linear viscous) is supported')
        if not(0 <= alpha <= 1): 
            raise ValueError('alpha should be between 0 and 1')

        ### Initialize
        mi   = np.zeros((3, 3, self.numdofs0)) # (xyz component, i-th vector, node)
        Eij  = np.zeros((6, self.numdofs0))    # (Eij component, node)
        lami = np.zeros((3, self.numdofs0))    # (i-th a^(2) eigenvalue, node)
        if len(ei_arg) == 3: 
            mi[:,:,:] = self.ei_tile(ei_arg) # set prescribed ei frame
        
        ### Calculate Eij etc. using specfabpy per node
        nlm = self.nlm_nodal()
        for nn in np.arange(self.numdofs0): 
            eigvecs, lami[:,nn] = sfcom.eigenframe(nlm[nn,:], symframe=self.symframe, modelplane=self.modelplane)
            if len(ei_arg) == 0:
                mi[:,0,nn], mi[:,1,nn], mi[:,2,nn] = eigvecs.T
            Eij[:,nn] = self.sf.Eij_tranisotropic(nlm[nn,:], mi[:,0,nn], mi[:,1,nn], mi[:,2,nn], Eij_grain,alpha,n_grain)
        
        """
        The enhancement-factor model depends on effective (homogenized) grain parameters, calibrated against deformation tests.
        For CPOs far from the calibration states, negative values *may* occur where Eij should tend to zero if truncation L is not large enough.
        """
        Eij[Eij < 0] = 1e-2 # Set negative E_ij to a very small value (flow inhibiting)
        
        ### Construct function spaces of quantities
        mi_df   = [Function(self.R30) for _ in range(3)] # (m1,m2,m3)
        Eij_df  = [Function(self.R0)  for _ in range(6)] # Eij enhancement tensor
        lami_df = [Function(self.R0)  for _ in range(3)] # a2 eigenvalues (lami)

        for ii in range(3): # m1,m2,m3
            lami_df[ii].vector()[:] = lami[ii,:]
            for jj in range(3): # x,y,z
                mi_df[ii].sub(jj).vector()[:] = mi[jj,ii,:]

        for kk in range(6): 
            Eij_df[kk].vector()[:] = Eij[kk,:]
            
        return mi_df, Eij_df, lami_df
        
    def ei_tile(self, ei):
        ei_tile  = np.zeros((3, 3, self.numdofs0)) # (xyz component, i-th vector, node)
        for ii in range(3): ei_tile[:,ii,:] = np.tile(ei[ii], (self.numdofs0,1)).T
        return ei_tile
        
