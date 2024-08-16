#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2024

"""
FEniCS interface for CPO dynamics using specfab
"""

import numpy as np
from dolfin import *
from ..specfabpy import specfabpy as sf__ # sf private copy 
from .. import common as sfcom

class CPO():
    
    """
    Class representing a single crystallographic axis
    """

    def __init__(self, mesh, boundaries, L, nu_multiplier=1, nu_realspace=1e-3, modelplane='xz', ds=None, nvec=None):

        ### Setup
        
        self.mesh, self.boundaries = mesh, boundaries
        self.ds = ds   if ds   is not None else Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        self.n  = nvec if nvec is not None else FacetNormal(self.mesh)
        self.L = int(L) # spectral truncation
        self.USE_REDUCED = True # use reduced representation of fabric state vector 

        self.modelplane = modelplane
        if self.modelplane not in ['xy', 'xz']: 
            raise ValueError('modelplane "%s" should be either "xy" or "xz"'%(self.modelplane))
                
        self.nu_realspace  = Constant(nu_realspace)  # Real-space stabilization (multiplicative constant of real-space Laplacian)
        self.nu_multiplier = Constant(nu_multiplier) # Multiplier of orientation-space regularization magnitude
                
        ### Initialize specfab
        
        self.sf = sf__
        self.lm, self.nlm_len_full = self.sf.init(self.L)
        self.nlm_len = self.sf.get_rnlm_len() if self.USE_REDUCED else self.nlm_len_full # use reduced or full form?
        
        ### Function spaces
        
        eletype, eleorder = 'CG', 1 # ordinary linear elements
        self.Rele = FiniteElement(eletype, self.mesh.ufl_cell(), eleorder)
        self.Vele = MixedElement([self.Rele]*self.nlm_len) # note: a VectorElement is nothing but a MixedElement combining (multiplying) "dim" copies of a FiniteElement
        self.R = FunctionSpace(self.mesh, self.Rele) 
        self.V = FunctionSpace(self.mesh, self.Vele)
        self.G = TensorFunctionSpace(self.mesh, eletype, eleorder, shape=(2,2)) # strain-rate and spin function space 
        self.numdofs = Function(self.R).vector().local_size() # for MPI support, get vector() size like this. Else could have used: numdofs = self.R.dim()

        if self.modelplane=='xy':
    
            self.Wele = MixedElement([self.Vele, self.Vele]) # real, imag components
            self.W = FunctionSpace(self.mesh, self.Wele)
            self.dofs_re = np.array([ self.W.sub(0).sub(ii).dofmap().dofs() for ii in range(self.nlm_len) ])
            self.dofs_im = np.array([ self.W.sub(1).sub(ii).dofmap().dofs() for ii in range(self.nlm_len) ])
        
            # Test and trial functions
            self.pr, self.pi = TrialFunctions(self.W) # unknowns (real, imag part of nlm coefs)
            self.qr, self.qi = TestFunctions(self.W)  # weight functions (real, imag part of nlm coefs)
            self.qr_sub, self.qi_sub = split(self.qr), split(self.qi) # for easy access of each subelement of the mixed element (real, imag parts)
        
        elif self.modelplane=='xz':    
                    
            self.Wele = self.Vele
            self.W = self.V
            self.dofs_re = np.array([ self.W.sub(ii).dofmap().dofs() for ii in range(self.nlm_len) ])
            
            # Test and trial functions
            self.pr = TrialFunction(self.W) # unknowns (real part of nlm coefs)
            self.qr = TestFunction(self.W)  # weight functions (real part of nlm coefs)
            self.pr_sub = split(self.pr)    # for easy access of each subelement
            self.qr_sub = split(self.qr)    # for easy access of each subelement
            
        ### Solution containers
        
        self.w      = Function(self.W) # Current solution
        self.w_prev = Function(self.W) # Previous solution
        
        ### Dynamical matrices
        
        self.nlm_dummy = np.zeros((self.nlm_len_full))
        self.Mk_LROT     = [ [Function(self.V) for ii in np.arange(self.nlm_len)] for ii in range(4) ] # rr, ri, ir, ii
        self.Mk_DDRX_src = [ [Function(self.V) for ii in np.arange(self.nlm_len)] for ii in range(4) ]
        self.Mk_CDRX     = [ [Function(self.V) for ii in np.arange(self.nlm_len)] for ii in range(4) ]
        self.Mk_REG      = [ [Function(self.V) for ii in np.arange(self.nlm_len)] for ii in range(4) ]
        
        ### Aux
        
        # Idealized states
        self.nlm_iso   = [1/np.sqrt(4*np.pi)] + [0]*(self.nlm_len-1) # Isotropic and normalized state
        self.nlm_zero  = [0]*(self.nlm_len)
                    

    def initialize(self, wr=None, wi=None):
        
        """
        Initialize uniform CPO field
        """
        
        if wr is None: wr, wi = self.nlm_iso, self.nlm_zero
        
        if self.modelplane=='xy':
            assign(self.w.sub(0), project(Constant(wr), self.V)) # real part
            assign(self.w.sub(1), project(Constant(wi), self.V)) # imag part
            
        elif self.modelplane=='xz':
            assign(self.w, project(Constant(wr), self.V)) # real part            


    def set_state(self, w, interp=False):
        if interp:
            raise ValueError('CPO.set_state() supports only setting function space vars, not interpolating expressions or constants.')
        else:
            self.w.assign(w)
            self.w_prev.assign(w)
                    

    def set_BCs(self, wr, wi, domid, domain=None):

        """
        Set boundary conditions
        """
    
        self.bcs = []
        if domain is None: domain = self.boundaries 
        
        for ii, did in enumerate(domid):
        
            if self.modelplane=='xy':
                self.bcs += [DirichletBC(self.W.sub(0), wr[ii], domain, did)] # real part
                self.bcs += [DirichletBC(self.W.sub(1), wi[ii], domain, did)] # imag part
                
            elif self.modelplane=='xz':
                self.bcs += [DirichletBC(self.W, wr[ii], domain, did)] # real part
                
    
    def set_isotropic_BCs(self, domid, domain=None):

        """
        Easy setting of isotropic boundary conditions
        """
    
        wr = [Constant(self.nlm_iso),]
        wi = [Constant(self.nlm_zero),]
        self.set_BCs(wr, wi, [domid,], domain=domain)
    

    def evolve(self, u, S, dt, iota=+1, Gamma0=None, Lambda0=None, steadystate=False):
    
        """
        Evolve CPO using Laplacian stabilized, Euler time integration
        """
    
        if self.w is None:
            raise ValueError('CPO state "w" not set. Did you forget to initialize the CPO field?')
            
        self.w_prev.assign(self.w) # current state (w) must be set
        F = self.weakform(u, S, dt, iota, Gamma0, Lambda0, steadystate=steadystate)
        solve(lhs(F)==rhs(F), self.w, self.bcs, solver_parameters={'linear_solver':'gmres', }) # fastest tested are: gmres, bicgstab, tfqmr (For non-symmetric problems, a Krylov solver for non-symmetric systems, such as GMRES, is a better choice)


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
        
        # Dynamical matrices at each DOF
        if ENABLE_LROT: M_LROT_nodal    = np.array([self.sf.reduce_M(self.sf.M_LROT(self.nlm_dummy, self.mat3d(Df[nn*4:(nn+1)*4]), self.mat3d(Wf[nn*4:(nn+1)*4]), iota, zeta), self.nlm_len) for nn in np.arange(self.numdofs)] )
        if ENABLE_DDRX: M_DDRX_src_nodal = np.array([self.sf.reduce_M(self.sf.M_DDRX_src(self.nlm_dummy, self.mat3d(Sf[nn*4:(nn+1)*4])), self.nlm_len) for nn in np.arange(self.numdofs)] )
        M_REG_nodal = np.array([self.sf.reduce_M(self.sf.M_REG(self.nlm_dummy, self.mat3d(Df[nn*4:(nn+1)*4])), self.nlm_len) for nn in np.arange(self.numdofs)] )
        
        # Populate entries of dynamical matrices 
        if   self.modelplane=='xy': krng = range(4) # rr, ri, ir, ii
        elif self.modelplane=='xz': krng = range(1) # rr
        for ii in np.arange(self.nlm_len):
            for kk in krng: 
                if ENABLE_LROT: self.Mk_LROT[kk][ii].vector()[:] = M_LROT_nodal[:,kk,ii,:].flatten()
                if ENABLE_DDRX: self.Mk_DDRX_src[kk][ii].vector()[:] = M_DDRX_src_nodal[:,kk,ii,:].flatten()
                self.Mk_REG[kk][ii].vector()[:] = M_REG_nodal[:,kk,ii,:].flatten()

        ### Construct weak form
        
        dtinv = Constant(1/dt)    
       
        if self.modelplane=='xy':

            # Real space advection, div(s*u)
            F  = dot(dot(u, nabla_grad(self.pr)), self.qr)*dx # real part
            F += dot(dot(u, nabla_grad(self.pi)), self.qi)*dx # imag part
            
            # Time derivative
            if not steadystate:
                F += dtinv * dot( (self.pr-self.w_prev.sub(0)), self.qr)*dx # real part
                F += dtinv * dot( (self.pi-self.w_prev.sub(1)), self.qi)*dx # imag part
           
            # Real space stabilization (Laplacian diffusion)
            F += self.nu_realspace * inner(grad(self.pr), grad(self.qr))*dx # real part
            F += self.nu_realspace * inner(grad(self.pi), grad(self.qi))*dx # imag part
     
            # Lattice rotation
            if ENABLE_LROT:
                Mrr_LROT, Mri_LROT, Mir_LROT, Mii_LROT = self.Mk_LROT # unpack for readability 
                F += -sum([ (dot(Mrr_LROT[ii], self.pr) + dot(Mri_LROT[ii], self.pi))*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # real part
                F += -sum([ (dot(Mir_LROT[ii], self.pr) + dot(Mii_LROT[ii], self.pi))*self.qi_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # imag part
            
            # Orientation space stabilization (hyper diffusion)
            Mrr_REG,  Mri_REG,  Mir_REG,  Mii_REG  = self.Mk_REG  # unpack for readability 
            F += -self.nu_multiplier * sum([ (dot(Mrr_REG[ii], self.pr) + dot(Mri_REG[ii], self.pi))*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # real part
            F += -self.nu_multiplier * sum([ (dot(Mir_REG[ii], self.pr) + dot(Mii_REG[ii], self.pi))*self.qi_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # imag part

        elif self.modelplane=='xz':
                
            # Real space advection, div(s*u)
            F = dot(dot(u, nabla_grad(self.pr)), self.qr)*dx # real part
                
            # Time derivative
            if not steadystate:
                F += dtinv * dot( (self.pr-self.w_prev), self.qr)*dx # real part
            
            # Real space stabilization (Laplacian diffusion)
            F += self.nu_realspace * inner(grad(self.pr), grad(self.qr))*dx # real part
     
            # Lattice rotation
            if ENABLE_LROT:
                Mrr_LROT, *_ = self.Mk_LROT # unpack for readability 
                F += -sum([ dot(Mrr_LROT[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # real part
            
            # DDRX 
            if ENABLE_DDRX:
                Mrr_DDRX_src, *_ = self.Mk_DDRX_src # unpack for readability 
                F_src  = sum([ -Gamma0*dot(Mrr_DDRX_src[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) 
                F_sink = sum([ -Gamma0*dot(Mrr_DDRX_src[0],self.w_prev)/self.w_prev[0]*self.pr_sub[ii] * self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # nonlinear sink term is linearized around previous solution (self.w_prev) following Rathmann and Lilien (2021)
#                F_sink = sum([ -Gamma0*dot(Mrr_DDRX_src[0],self.pr)*self.pr_sub[ii] * self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # nonlinear sink term is linearized around previous solution (self.w_prev) following Rathmann and Lilien (2021)
                F += F_src - F_sink
            
            # Orientation space stabilization (hyper diffusion)
            Mrr_REG,  *_ = self.Mk_REG  # unpack for readability 
            F += -self.nu_multiplier * sum([ dot(Mrr_REG[ii], self.pr)*self.qr_sub[ii]*dx for ii in np.arange(self.nlm_len)]) # real part

        return F
        
        
    def mat3d(self, mat2d): 
        return sfcom.mat3d(mat2d, self.modelplane, reshape=True) # from common.py
       
        
    def get_nlm(self, x, y):
    
        """
        Extract CPO state gridpoint-wise (full form of nlm)
        """
        
        if   self.modelplane=='xy': nlm = self.w.sub(0)(x,y) + 1j * self.w.sub(1)(x,y)
        elif self.modelplane=='xz': nlm = self.w(x,y) + 0j
        return self.sf.rnlm_to_nlm(nlm, self.nlm_len_full) if self.USE_REDUCED else nlm
        
       
    def eigenframe(self, x, y, modelplane=None, extrapolate=False):

        """
        Extract CPO eigenvectors and eigenvalues for selected locations (x,y)
        """
        
        self.w.set_allow_extrapolation(extrapolate)

        xf, yf = np.array(x, ndmin=1), np.array(y, ndmin=1)
        if len(xf.shape) > 1: raise ValueError('eigenframe(): only 1D (x,y) are supported (i.e. flattened)')
        
        N = len(xf)
        eigvecs, eigvals = np.zeros((N,3,3)), np.zeros((N,3))
        for ii in np.arange(N):
            a2 = self.sf.a2(self.get_nlm(xf[ii],yf[ii]))
            eigvecs[ii,:,:], eigvals[ii,:] = sfcom.eigenframe(a2, modelplane=modelplane)

        return (eigvecs[0,:,:], eigvals[0,:]) if N==1 else (eigvecs, eigvals)
                   
        
    def apply_bounds(self):
    
        """ 
        Renormalize power spectrum if it exceeds that of the delta function (numerical overshoot)
        """
        
        w  = self.w.vector()[:]    # unbounded solution
        wb = np.zeros(np.shape(w)) # bounded solution
        
        for nn in np.arange(self.numdofs): 

            if   self.modelplane=='xy': nlm = w[self.dofs_re[:,nn]] + 1j*w[self.dofs_im[:,nn]]
            elif self.modelplane=='xz': nlm = w[self.dofs_re[:,nn]] + 0j
                
            if self.USE_REDUCED:
                nlm_full = self.sf.rnlm_to_nlm(nlm, self.nlm_len_full)
                nlm_full_bnd = self.sf.apply_bounds(nlm_full)
                nlm_bnd = self.sf.nlm_to_rnlm(nlm_full_bnd, self.nlm_len)
            else:
                nlm_bnd = self.sf.apply_bounds(nlm)
                
            if self.modelplane=='xy': 
                wb[self.dofs_re[:,nn]] = np.real(nlm_bnd)
                wb[self.dofs_im[:,nn]] = np.imag(nlm_bnd)  
                
            elif self.modelplane=='xz': 
                wb[self.dofs_re[:,nn]] = np.real(nlm_bnd)
            
        self.w.vector()[:] = wb.copy()
                          
