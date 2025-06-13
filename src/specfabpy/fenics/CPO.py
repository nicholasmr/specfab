#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>

"""
FEniCS interface for CPO dynamics using specfab
"""

import code, copy # code.interact(local=locals())
import numpy as np
from datetime import datetime
from dolfin import *
from ..specfabpy import specfabpy as sf__ # sf private copy 
from .. import common as sfcom

class CPO():
    
    """
    Class representing a single crystallographic axis
    """

    def __init__(self, mesh, boundaries, L, nu_multiplier=1, nu_realspace=1e-3, modelplane='xz', symframe=-1, ds=None, nvec=None):

        ### Check args
        
        if modelplane not in ['xz', 'xy']:
            raise ValueError('modelplane "%s" not supported, must be "xz"'%(modelplane))

        ### Setup
        
        self.mesh, self.boundaries = mesh, boundaries
        self.ds = ds   if ds   is not None else Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        self.n  = nvec if nvec is not None else FacetNormal(self.mesh)
        self.L = int(L) # spectral truncation
        self.USE_REDUCED = True # use reduced representation of fabric state vector 
        self.symframe = symframe
        self.modelplane = modelplane
        self.nu_realspace  = nu_realspace  # Real-space stabilization (multiplicative constant of real-space Laplacian)
        self.nu_multiplier = nu_multiplier # Multiplier of orientation-space regularization magnitude
                
        ### Initialize fortran module
        
        self.sf = sf__
        self.lm, self.nlm_len_full = self.sf.init(self.L)
        self.nlm_len = self.sf.get_rnlm_len() if self.USE_REDUCED else self.nlm_len_full # use reduced or full form?
        
        ### Function spaces for CPO dynamics
        
        ele = ('CG',1) # ordinary linear elements
        eletype, eleorder = ele
        self.Rele = FiniteElement(eletype, self.mesh.ufl_cell(), eleorder)
        self.Vele = MixedElement([self.Rele]*self.nlm_len) # note: a VectorElement is nothing but a MixedElement combining (multiplying) "dim" copies of a FiniteElement
        self.R = FunctionSpace(self.mesh, self.Rele) 
        self.V = FunctionSpace(self.mesh, self.Vele)
        self.G = TensorFunctionSpace(self.mesh, *ele, shape=(2,2)) # strain-rate and spin function space 
        self.numdofs = Function(self.R).vector().local_size() # for MPI support, get vector() size like this. Else could have used: numdofs = self.R.dim()

        # Ranges
        self.dofs  = range(self.numdofs)
        self.srng  = range(self.nlm_len)

        if self.modelplane=='xy':
    
            self.Sele = MixedElement([self.Vele, self.Vele]) # real, imag components
            self.S = FunctionSpace(self.mesh, self.Sele)
            self.dofs_re = np.array([ self.S.sub(0).sub(ii).dofmap().dofs() for ii in self.srng ])
            self.dofs_im = np.array([ self.S.sub(1).sub(ii).dofmap().dofs() for ii in self.srng ])
        
            # Test and trial functions
            self.sr, self.si = TrialFunctions(self.S) # unknowns (real, imag part of nlm coefs)
            self.wr, self.wi = TestFunctions(self.S)  # weight functions (real, imag part of nlm coefs)
            self.wr_sub, self.wi_sub = split(self.wr), split(self.wi) # for easy access of each subelement of the mixed element (real, imag parts)
        
        elif self.modelplane=='xz':    
                    
            self.Sele = self.Vele
            self.S = self.V
            self.dofs_re = np.array([ self.S.sub(ii).dofmap().dofs() for ii in self.srng ])
            
            # Test and trial functions
            self.sr = TrialFunction(self.S) # unknowns (real part of nlm coefs)
            self.wr = TestFunction(self.S)  # weight functions (real part of nlm coefs)
            self.sr_sub = split(self.sr)    # for easy access of each subelement
            self.wr_sub = split(self.wr)    # for easy access of each subelement

        ### Solution containers

        self.s  = Function(self.S) # Current solution
        self.s0 = Function(self.S) # Previous solution

        ### Dynamical matrices

        self.nlm_dummy = np.zeros((self.nlm_len_full))
        self.Mk_LROT     = [ [Function(self.V) for ii in self.srng] for ii in range(4) ] # rr, ri, ir, ii
        self.Mk_DDRX_src = [ [Function(self.V) for ii in self.srng] for ii in range(4) ]
        self.Mk_CDRX     = [ [Function(self.V) for ii in self.srng] for ii in range(4) ]
        self.Mk_REG      = [ [Function(self.V) for ii in self.srng] for ii in range(4) ]

        ### Aux

        # Idealized states
        self.nlm_iso   = [1/np.sqrt(4*np.pi)] + [0]*(self.nlm_len-1) # Isotropic and normalized state
        self.nlm_zero  = [0]*(self.nlm_len)
        self.M_zero    = np.zeros((self.numdofs, self.nlm_len, self.nlm_len))
        

    def initialize(self, sr=None, si=None):
        """
        Initialize uniform CPO field
        """
        if sr is None: 
            sr, si = self.nlm_iso, self.nlm_zero
        s = Function(self.S)
        if self.modelplane=='xy':
            assign(s.sub(0), project(Constant(sr), self.V)) # real part
            assign(s.sub(1), project(Constant(si), self.V)) # imag part
        elif self.modelplane=='xz':
            assign(s, project(Constant(sr), self.V)) # real part
        self.set_state(s)

    def set_state(self, s, interp=False):
        if interp:
            raise ValueError('CPO.set_state() supports only setting function space vars, not interpolating expressions or constants.')
        else:
            self.s.assign(s)
            self.s0.assign(s)

    def set_BCs(self, sr, si, domids, domain=None):
        """
        Set boundary conditions
        """    
        self.bcs = []
        if domain is None: domain = self.boundaries 
        for ii, did in enumerate(domids):
            if self.modelplane=='xy':
                self.bcs += [DirichletBC(self.S.sub(0), sr[ii], domain, did)] # real part
                self.bcs += [DirichletBC(self.S.sub(1), si[ii], domain, did)] # imag part
            elif self.modelplane=='xz':
                self.bcs += [DirichletBC(self.S, sr[ii], domain, did)] # real part
                
    def set_isotropic_BCs(self, domids, domain=None):
        """
        Easy setting of isotropic boundary conditions
        """
        sr = [Constant(self.nlm_iso)]  * len(domids)
        si = [Constant(self.nlm_zero)] * len(domids)
        self.set_BCs(sr, si, domids, domain=domain)

    def evolve(self, u, S, dt, iota=+1, Gamma0=None, Lambda0=None, disable_advection=False):
        """
        The fabric solver, called to step fabric field forward in time by amount dt
        """    
        self.s0.assign(self.s) # current state self.s must be set
        s = (self.sr,None) if self.modelplane=='xz' else (self.sr,self.si) # (sr, si)
        F = self._weakform(*s, u,S, dt, iota,Gamma0,Lambda0, disable_advection=disable_advection)
        solve(lhs(F)==rhs(F), self.s, self.bcs, solver_parameters={'linear_solver':'gmres'}) # fastest tested are: gmres, bicgstab, tfqmr --- note this is a non-symmetric system!

    def solvesteady(self, u, S, iota=+1, Gamma0=None, Lambda0=None, LROT_guess=False, nu_realspacemul=None, **kwargs):

        dt = None # select steady weak form

        # LROT solution (+CDRX+REG) *or* used as init guess for nonlinear problem
        if Gamma0 is None or LROT_guess:
            print('*** Solving linear LROT problem %s'%('(used as initial guess)' if (LROT_guess and (Gamma0 is not None)) else ''))
            self.evolve(u,S, dt, iota=iota, Gamma0=None, Lambda0=Lambda0, **kwargs) # solution stored in self.s 
            
        # LROT + DDRX solution (+CDRX+REG)
        if Gamma0 is not None: # nonlinear problem?
            print('*** Solving nonlinear LROT+DDRX problem')
            if not isinstance(Gamma0, list): Gamma0 = [Gamma0,]
            nu_realspace0 = copy.deepcopy(self.nu_realspace)
            for ii, gam0 in enumerate(Gamma0):
                self.nu_realspace = nu_realspacemul[ii]*nu_realspace0 # gradually refine regularization, too
                print('...approaching solution w/ nu_realmul=%.1f (step %i of %i)'%(nu_realspacemul[ii], ii+1, len(Gamma0)))
                s = (self.s, None) if self.modelplane=='xz' else (self.s.sub(0), self.s.sub(1)) # (sr, si)
                F = self._weakform(*s, u,S, dt, iota,gam0,Lambda0, **kwargs)
                solve(F==0, self.s, self.bcs, solver_parameters={'newton_solver':{'linear_solver':'gmres', 'preconditioner':'none'}})

    def _weakform(self, sr,si, u,S, dt, iota, Gamma0, Lambda0, zeta=0, disable_advection=False):
        
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
        
        # Dynamical matrices M_* at each node as np arrays (indexing is M[node,kk,row,column])
        M_LROT     = np.array([self._M_reduced(self.sf.M_LROT,     self.nlm_dummy, D3[nn], W3[nn], iota, zeta) for nn in self.dofs] ) if ENABLE_LROT else self.M_zero
        M_DDRX_src = np.array([self._M_reduced(self.sf.M_DDRX_src, self.nlm_dummy, S3[nn])                     for nn in self.dofs] ) if ENABLE_DDRX else self.M_zero
        M_REG      = np.array([self._M_reduced(self.sf.M_REG,      self.nlm_dummy, D3[nn])                     for nn in self.dofs] ) if ENABLE_REG  else self.M_zero
        
        # Populate entries of dynamical matrices 
        if   self.modelplane=='xy': krng = range(4) # rr, ri, ir, ii
        elif self.modelplane=='xz': krng = range(1) # rr
        for ii in self.srng:
            for kk in krng: 
                if ENABLE_LROT: self.Mk_LROT[kk][ii].vector()[:]     = M_LROT[:,kk,ii,:].flatten()
                if ENABLE_DDRX: self.Mk_DDRX_src[kk][ii].vector()[:] = M_DDRX_src[:,kk,ii,:].flatten()
                if ENABLE_REG:  self.Mk_REG[kk][ii].vector()[:]      = M_REG[:,kk,ii,:].flatten()

        ### Construct weak form
        
        if self.modelplane=='xy':
            
            # Real space stabilization (Laplacian diffusion)
            F  = Constant(self.nu_realspace) * inner(grad(sr), grad(self.wr))*dx # real part
            F += Constant(self.nu_realspace) * inner(grad(si), grad(self.wi))*dx # imag part

            # Time derivative
            if dt is not None:
                F += Constant(1/dt) * dot( (sr-self.s0.sub(0)), self.wr)*dx # real part
                F += Constant(1/dt) * dot( (si-self.s0.sub(1)), self.wi)*dx # imag part

            # Real space advection
            if not disable_advection:
                F += dot(dot(u, nabla_grad(sr)), self.wr)*dx # real part
                F += dot(dot(u, nabla_grad(si)), self.wi)*dx # imag part
     
            # Lattice rotation
            if ENABLE_LROT:
                Mrr_LROT, Mri_LROT, Mir_LROT, Mii_LROT = self.Mk_LROT # unpack for readability 
                F += -sum([ (dot(Mrr_LROT[ii], sr) + dot(Mri_LROT[ii], si))*self.wr_sub[ii]*dx for ii in self.srng]) # real part
                F += -sum([ (dot(Mir_LROT[ii], sr) + dot(Mii_LROT[ii], si))*self.wi_sub[ii]*dx for ii in self.srng]) # imag part
            
            # DDRX 
            if ENABLE_DDRX:
                raise ValueError('CPO(): DDRX evolution is not yet supported for modelplane="xy"')
            
            # Orientation space stabilization (hyper diffusion)
            if ENABLE_REG:
                Mrr_REG,  Mri_REG,  Mir_REG,  Mii_REG  = self.Mk_REG  # unpack for readability 
                F += -Constant(self.nu_multiplier) * sum([ (dot(Mrr_REG[ii], sr) + dot(Mri_REG[ii], si))*self.wr_sub[ii]*dx for ii in self.srng]) # real part
                F += -Constant(self.nu_multiplier) * sum([ (dot(Mir_REG[ii], sr) + dot(Mii_REG[ii], si))*self.wi_sub[ii]*dx for ii in self.srng]) # imag part

        elif self.modelplane=='xz':
                
            # dummy zero term to make rhs(F) work when solving steady state problem 
            # this can probably be removed once the SSA source/sink terms are added
            s_null = Function(self.S)
            s_null.vector()[:] = 0.0
            F = dot(s_null, self.wr)*dx

            # Time derivative
            if dt is not None:
                F += Constant(1/dt) * dot( (sr-self.s0), self.wr)*dx # real part

            # Real space advection
            if not disable_advection:
                F += dot(dot(u, nabla_grad(sr)), self.wr)*dx # real part
                
            # Real space stabilization (Laplacian diffusion)
            F += Constant(self.nu_realspace) * inner(grad(sr), grad(self.wr))*dx # real part
    
            # Lattice rotation
            if ENABLE_LROT:
                Mrr_LROT, *_ = self.Mk_LROT # unpack for readability 
                F += -sum([ dot(Mrr_LROT[ii], sr)*self.wr_sub[ii]*dx for ii in self.srng]) # real part
            
            # DDRX 
            if ENABLE_DDRX:
                Mrr_DDRX_src, *_ = self.Mk_DDRX_src # unpack for readability 
                F_src = -sum([Gamma0*dot(Mrr_DDRX_src[ii], sr)*self.wr_sub[ii]*dx for ii in self.srng]) 
                
                sr_sub = split(sr) if dt is None else self.sr_sub
                if dt is None: 
                    # Steady state solver? Use nonlinear iteration for DDRX-activated problem
                    D = dot(Mrr_DDRX_src[0], sr)/sr_sub[0] # <D> (Rathmann et al., 2025)
                else: 
                    # Time-dependent probem? Linearized problem by linearizing <D>
                    D = dot(Mrr_DDRX_src[0], self.s0)/self.s0[0] # <D> estimated using previous solution (Rathmann and Lilien, 2021)
                F_snk = -sum([Gamma0*D*sr_sub[ii]*self.wr_sub[ii]*dx for ii in self.srng])
                
                F += F_src - F_snk
            
            # Orientation space stabilization (hyper diffusion)
            if ENABLE_REG:
                Mrr_REG,  *_ = self.Mk_REG  # unpack for readability 
                F += -Constant(self.nu_multiplier) * sum([ dot(Mrr_REG[ii], sr)*self.wr_sub[ii]*dx for ii in self.srng]) # real part

        return F
        
    def _M_reduced(self, Mfun, *args, **kwargs):
        return self.sf.reduce_M(Mfun(*args, **kwargs), self.nlm_len) # full form -> reduced form of M_*

    def mat3d(self, D2, dn=4):
        D2 = np.array([ D2[nn*dn:(nn+1)*dn] for nn in self.dofs ]) # to [node, D2 flat index 0 to 3]
        return sfcom.mat3d_arr(D2, self.modelplane, reshape=True) # [node,3,3]
        
    def eigenframe(self, *p, reduce2D=False, **kwargs):
        """
        Extract a2 eigenvectors (mi) and eigenvalues (lami) at specified list of points p=([x1,x2,...], [y1,y2,...])
        Note that you can back-construct a2 = lam1*np.outer(m1,m1) + lam2*np.outer(m2,m2) + lam3*np.outer(m3,m3) 
        """
        xf, yf = np.array(p[0], ndmin=1), np.array(p[1], ndmin=1) # (point num, x/y value)
        nlm = np.array([self.get_nlm(xf[pp],yf[pp], **kwargs) for pp in range(len(xf))]) # (point num, nlm coefs)
        mi, lami = sfcom.eigenframe(nlm, symframe=self.symframe, modelplane=self.modelplane)
        
        if reduce2D:
            if nlm.shape[0] > 1:
                if self.modelplane == 'xz': mi, lami = mi[:,:,[0,2]], lami[:,[0,2]]
                if self.modelplane == 'xy': mi, lami = mi[:,:,[0,1]], lami[:,[0,1]]
            else:
                if self.modelplane == 'xz': mi, lami = mi[:,[0,2]], lami[[0,2]]
                if self.modelplane == 'xy': mi, lami = mi[:,[0,1]], lami[[0,1]]
                
        return (mi, lami)

    def get_nlm(self, *p, xz2xy=False):    
        """
        Extract CPO state gridpoint-wise (full form of nlm)
        """
        if   self.modelplane=='xy': nlm_ = self.s.sub(0)(*p) + 1j * self.s.sub(1)(*p)
        elif self.modelplane=='xz': nlm_ = self.s(*p) + 0j
        nlm = self.sf.rnlm_to_nlm(nlm_, self.nlm_len_full) if self.USE_REDUCED else nlm_
        return self.sf.rotate_nlm_xz2xy(nlm) if xz2xy else nlm 
                
    def apply_bounds(self):
        """ 
        Renormalize power spectrum if it exceeds that of the delta function (numerical overshoot)
        """
        
        s  = self.s.vector()[:]    # unbounded solution
        sb = np.zeros(np.shape(s)) # bounded solution
        
        for nn in self.dofs: 

            if   self.modelplane=='xy': nlm = s[self.dofs_re[:,nn]] + 1j*s[self.dofs_im[:,nn]]
            elif self.modelplane=='xz': nlm = s[self.dofs_re[:,nn]] + 0j
                
            if self.USE_REDUCED:
                nlm_full = self.sf.rnlm_to_nlm(nlm, self.nlm_len_full)
                nlm_full_bnd = self.sf.apply_bounds(nlm_full)
                nlm_bnd = self.sf.nlm_to_rnlm(nlm_full_bnd, self.nlm_len)
            else:
                nlm_bnd = self.sf.apply_bounds(nlm)
                
            if self.modelplane=='xy': 
                sb[self.dofs_re[:,nn]] = np.real(nlm_bnd)
                sb[self.dofs_im[:,nn]] = np.imag(nlm_bnd)  
            elif self.modelplane=='xz': 
                sb[self.dofs_re[:,nn]] = np.real(nlm_bnd)
            
        self.s.vector()[:] = sb.copy()
                          
