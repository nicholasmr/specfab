#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2024

"""
FEniCS interface for calculating bulk enhancement factors given CPO field
"""

import numpy as np
from dolfin import *
from ..specfabpy import specfabpy as sf__ # sf private copy 
from ..common import *
from .rheology import Orthotropic, Isotropic

import copy, sys, time, code # code.interact(local=locals())

class EnhancementFactor():

    def __init__(self, mesh, L, modelplane='xz'):

        ### Setup

        self.mesh = mesh
        self.L = L        
        self.modelplane = modelplane
        self.USE_REDUCED = True # use reduced representation of fabric state vector 
        
        self.sf = sf__
        self.lm, self.nlm_len_full = self.sf.init(self.L)
        self.nlm_len = self.sf.get_rnlm_len() if self.USE_REDUCED else self.nlm_len_full # use reduced or full form?
        
        ### Function spaces
        
        # Store CPO information element-wise: node located at element center, CPO constant over element
        eletype, eleorder = 'DG', 0 
        self.Rele = FiniteElement(eletype, self.mesh.ufl_cell(), eleorder)
        self.R  = FunctionSpace(self.mesh, self.Rele) # for scalars 
        self.R3 = VectorFunctionSpace(self.mesh, eletype, eleorder, dim=3) # for vectors
        self.numdofs = Function(self.R).vector().local_size()
        self.R3_assigner = FunctionAssigner(self.R3, [self.R, self.R, self.R]) # for constructing vectors from their components
        self.G = TensorFunctionSpace(self.mesh, eletype, eleorder, shape=(2,2)) # strain-rate and spin function space 
        

    def V2R(self, wV):
        wV_sub = wV.split()
        wR_sub = [project(wV_sub[ii], self.R) for ii in range(self.nlm_len)] 
        return wR_sub
        
        
    def wR(self, w):
        if   self.modelplane=='xy': return self.V2R(w.sub(0)), self.V2R(w.sub(1)) # wr_R, wi_R
        elif self.modelplane=='xz': return self.V2R(w), [project(Constant(0), self.R)]*self.nlm_len # wr_R, 0
        
       
    def nlm_nodal(self, w):
        wr_R, wi_R = self.wR(w)
        rnlm = np.array([ wr_R[ii].vector()[:] + 1j*wi_R[ii].vector()[:] for ii in range(self.nlm_len) ])
        nlm = np.array([ self.sf.rnlm_to_nlm(rnlm[:,nn], self.nlm_len_full) for nn in range(self.numdofs) ])
        return nlm.T # nlm[coef,node]

 
    def get_nlm(self, w, x,y):
        wr_R, wi_R = self.wR(w)
        rnlm = np.array([ wr_R[ii].vector()[:] + 1j*wi_R[ii].vector()[:] for ii in range(self.nlm_len) ])
        return self.sf.rnlm_to_nlm(rnlm, self.nlm_len_full)
        
        
    def ei_tile(self, ei):
        ei_tile  = np.zeros((3, 3, self.numdofs)) # (xyz component, i-th vector, node)
        for ii in range(3): ei_tile[:,ii,:] = np.tile(ei[ii], (self.numdofs,1)).T
        return ei_tile            
        
        
    def np_to_func(self, ei, Eij, ai):

        # CPO eigenvectors
        ei_df  = [Function(self.R3) for _ in range(3)] # list of ei
        eij_df = [Function(self.R)  for _ in range(3)] # list of z,y,z components of a given ei
        for ii in range(3): # loop over vectors (e1,e2,e3)
            for jj in range(3): # loop over vector components
                eij_df[jj].vector()[:] = ei[jj,ii,:]
            self.R3_assigner.assign(ei_df[ii], eij_df) # set vector field components
        
        # Enhancement factors
        Eij_df = [Function(self.R) for _ in range(6)]
        for kk in range(6): Eij_df[kk].vector()[:] = Eij[kk,:]
        
        # Eigenvalues
        ai_df = [Function(self.R) for _ in range(3)]
        for ii in range(3): ai_df[ii].vector()[:] = ai[ii,:]
        
        return ei_df, Eij_df, ai_df
        
        
    def Eij_tranisotropic(self, w, Eij_grain, alpha, n_grain, ei_arg=()):
   
        """
        Bulk enhancement factors w.r.t. ei=(e1,e2,e3) axes for *transversely isotropic* grains.
        If ei=() then a^(2) eigenvectors are taken, ei=(m1,m2,m3), and Eij are the eigenenhancements.
        """
    
        if n_grain != 1: raise ValueError('only n_grain = 1 (linear viscous) is supported')
        if not(0 <= alpha <= 1): raise ValueError('alpha should be between 0 and 1')
       
        ei  = np.zeros((3, 3, self.numdofs)) # (xyz component, i-th vector, node)
        Eij = np.zeros((6, self.numdofs))    # (Eij component, node)
        ai  = np.zeros((3, self.numdofs))    # (i-th a^(2) eigenvalue, node)

        if len(ei_arg) == 3: ei[:,:,:] = self.ei_tile(ei_arg) # set prescribed ei frame
        nlm = self.nlm_nodal(w) # (nlm component, node)
        
        for nn in np.arange(self.numdofs): 

            eigvecs, eigvals = eigenframe(self.sf.a2(nlm[:,nn]), modelplane=self.modelplane)
            ai[:,nn] = eigvals
            
            if len(ei_arg) == 0:
                ei[:,0,nn], ei[:,1,nn], ei[:,2,nn] = eigvecs.T # principal directions of a^(2)
                
            Eij[:,nn] = self.sf.Eij_tranisotropic(nlm[:,nn], ei[:,0,nn], ei[:,1,nn], ei[:,2,nn], Eij_grain,alpha,n_grain) # 3x3 tensor, eigenenhancements
        
        """
        The enhancement-factor model depends on effective (homogenized) grain parameters, calibrated against deformation tests.
        For CPOs far from the calibration states, negative values may occur where Eij should tend to zero if truncation L is not large enough.
        """
        Eij[Eij < 0] = 1e-2 # Set negative E_ij to a very small value (flow inhibiting)
        
        return self.np_to_func(ei, Eij, ai) # return ei_df, Eij_df, ai_df
        
    
    def Eij_orthotropic(self, wb, wn, Eij_grain, alpha, n_grain, ei_arg=()):
   
        """
        Bulk enhancement factors w.r.t. ei=(e1,e2,e3) axes for *orthotropic* grains.
        If ei=() then a^(2) eigenvectors are taken, ei=(m1,m2,m3), and Eij are the eigenenhancements.
        """
    
        if n_grain != 1: raise ValueError('only n_grain = 1 (linear viscous) is supported')
        if not(0 <= alpha <= 1): raise ValueError('alpha should be between 0 and 1')
       
        ei  = np.zeros((3, 3, self.numdofs)) # (xyz component, i-th vector, node)
        Eij = np.zeros((6, self.numdofs))    # (Eij component, node)
        ai  = np.zeros((3, self.numdofs))    # (i-th a^(2) eigenvalue, node)

        if len(ei_arg) == 3: ei[:,:,:] = self.ei_tile(ei_arg) # set prescribed ei frame
        blm = self.nlm_nodal(wb) # (blm component, node)
        nlm = self.nlm_nodal(wn) # (nlm component, node)
        vlm = 0*nlm[:,0] # calculate from (blm,nlm) using joint ODF
        
        for nn in np.arange(self.numdofs): 

            eigvecs, eigvals = eigenframe(self.sf.a2(nlm[:,nn]), modelplane=self.modelplane)
            ai[:,nn] = eigvals
            
            if len(ei_arg) == 0:
                ei[:,0,nn], ei[:,1,nn], ei[:,2,nn] = eigvecs.T # principal directions of a^(2)
                
            Eij[:,nn] = self.sf.Eij_orthotropic(blm[:,nn], nlm[:,nn], vlm, ei[:,0,nn], ei[:,1,nn], ei[:,2,nn], Eij_grain,alpha,n_grain) # 3x3 tensor, eigenenhancements
        
            if np.any(Eij[:,nn] < 0): 
                print('[!!] Correction needed @ %03i :: Eij = '%(nn), Eij[:,nn])
#                print('... m1 = ', ei[:,0,nn])
#                print('... m2 = ', ei[:,1,nn])
                print('... m3 = ', ei[:,2,nn])
        
        """
        The enhancement-factor model depends on effective (homogenized) grain parameters, calibrated against deformation tests.
        For CPOs far from the calibration states, negative values may occur where Eij should tend to zero if truncation L is not large enough.
        """
        Eij[Eij < 0] = 1e-2 # Set negative E_ij to a very small value (flow inhibiting)
        
        return self.np_to_func(ei, Eij, ai) # return ei_df, Eij_df, ai_df
        
    
    def E_EIE(self, u, Eij, mi, n, q=None, dim=2):
    
        """
        Equivalent isotropic enhancement factor (EIE) estimated from tensorial viscous structure.
        """
        
        if q is None: q = -n-1 # ideal solution
        
        D = sym(grad(u))
        if dim==2: D = as_tensor(mat3d(D, self.modelplane)) # 2x2 to 3x3 strain-rate tensor
        
        ort = Orthotropic(n=Constant(n))
        iso = Isotropic(  n=Constant(n))

        # Effective strain rates
        epsE_ort = sqrt(inner(ort.C_inv(D, mi, Eij), D))
        epsE_iso = sqrt(inner(iso.C_inv(D),D))

        # Equivalent isotropic enhancement factor
        E = (epsE_ort/epsE_iso)**(q)
        
        return project(E, self.R) # E as function
        

    def E_CAFFE(self, u, w, Emin, Emax):
    
        """
        CAFFE model (Placidi et al., 2010)

        Assumes D(stress tensor) = D(strain-rate tensor) where D is the deformability
        """

        Df = project( sym(grad(u)), self.G).vector()[:] # flattened strain-rate tensor
        E = np.zeros((self.numdofs))
        nlm = self.nlm_nodal(w) # (nlm component, node)
        
        for nn in np.arange(self.numdofs): 
            D = mat3d(Df[nn*4:(nn+1)*4], self.modelplane, reshape=True) # strain-rate tensor of nn-th node
            E[nn] = sf__.E_CAFFE(D, nlm[:,nn], Emin, Emax)
        
        E_df = Function(self.R)
        E_df.vector()[:] = E[:]
        
        return E_df


