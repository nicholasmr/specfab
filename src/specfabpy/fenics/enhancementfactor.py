#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>

"""
FEniCS interface for calculating bulk enhancement factors given a CPO field
"""

import numpy as np
#import code # code.interact(local=locals())
from dolfin import *
from ..specfabpy import specfabpy as sf__ # sf private copy 
from .. import common as sfcom
from .rheology import Orthotropic, Isotropic

class EnhancementFactor():

    def __init__(self, mesh, L, enhancementmodel='LTS', homoparams=(), modelplane='xz', symframe=-1, ele=()):

        ### Setup

        self.mesh = mesh
        self.L = L
        self.enhancementmodel = enhancementmodel # 'LTS', 'APEX', ...
        self.homoparams = homoparams
        self.modelplane = modelplane
        self.symframe = symframe # rheological symmetry frame: -1 is a2, 0-5 are the a4 eigentensors
        self.USE_REDUCED = True # use reduced representation of fabric state vector 
        
        if self.enhancementmodel not in ['LTS', 'APEX']:           raise ValueError('Invalid enhancementmodel "%s"'%(self.enhancementmodel))
        if self.enhancementmodel=='LTS'  and len(homoparams) != 3: raise ValueError('LTS model requires homoparams = (Eij_grain, alpha, n_grain)')
        if self.enhancementmodel=='APEX' and len(homoparams) != 2: raise ValueError('APEX model requires homoparams = (Emin, Emax)')
        
        self.sf = sf__
        self.lm, self.nlm_len_full = self.sf.init(self.L)
        self.nlm_len = self.sf.get_rnlm_len() if self.USE_REDUCED else self.nlm_len_full # use reduced or full form?
        
        ### FE spaces for viscous anisotropy, a2 eigenvalues, J index, ...
        
        if len(ele)==0: ele = ('DG',0) # WARNING: DG0 is the only tested safe element type for forward modeling of flow problems

        self.R = FunctionSpace(self.mesh, *ele)
        self.V = VectorFunctionSpace(self.mesh, *ele, dim=3) # for vectors
        self.V_assigner = FunctionAssigner(self.V, [self.R, self.R, self.R]) # for constructing vectors from their components
        self.G = TensorFunctionSpace(self.mesh, *ele, shape=(2,2))
        self.T = TensorFunctionSpace(self.mesh, *ele, shape=(3,3))
        self.numdofs = self.R.dim()
        
        ### Aux/shortcuts
        
        self.dofs0 = np.arange(self.numdofs)
        self.srng  = np.arange(self.nlm_len)


    def get_nlm(self, s, x,y):
        """
        Get numpy nlm array at point (x,y) for state vector field s
        """
        sr_R, si_R = self._sR(s)
        rnlm = np.array([ sr_R[ii].vector()[:] + 1j*si_R[ii].vector()[:] for ii in self.srng ])
        return self.sf.rnlm_to_nlm(rnlm, self.nlm_len_full)

    def _nlm_nodal(self, s):
        sr_R, si_R = self._sR(s)
        rnlm = np.array([ sr_R[ii].vector()[:] + 1j*si_R[ii].vector()[:] for ii in self.srng ])
        nlm  = np.array([ self.sf.rnlm_to_nlm(rnlm[:,nn], self.nlm_len_full) for nn in self.dofs0 ])
        return nlm # nlm[node,component]

    def _sR(self, s):
        if   self.modelplane=='xy': return self._V2R(s.sub(0)), self._V2R(s.sub(1)) # sr_R, si_R
        elif self.modelplane=='xz': return self._V2R(s), [project(Constant(0), self.R)]*self.nlm_len # sr_R, 0
                
    def _V2R(self, sV):
        sV_sub = sV.split()
        sR_sub = [project(sV_sub[ii], self.R) for ii in self.srng] 
        return sR_sub

    def _np2func(self, ei, Eij, ai):
    
        mi_df   = [Function(self.V) for _ in range(3)] # (m1,m2,m3) fabric principal directions
        Eij_df  = [Function(self.R) for _ in range(6)] # Eij enhancement tensor
        lami_df = [Function(self.R) for _ in range(3)] # a2 eigenvalues (lami)
        eij_df  = [Function(self.R) for _ in range(3)] # j-th component of ei
                
        for ii in range(3): # loop over vectors (e1,e2,e3)
            for jj in range(3): # loop over vector components of ei
                eij_df[jj].vector()[:] = ei[ii][:,jj] # [i,node,xyz]
            self.V_assigner.assign(mi_df[ii], eij_df) # set vector field components
            
        for kk in range(6): 
            Eij_df[kk].vector()[:] = Eij[:,kk]
        
        for ii in range(3): 
            lami_df[ii].vector()[:] = ai[:,ii] 
        
        return (mi_df, Eij_df, lami_df)
        
    def mat3d(self, D2, dn=4):
        D2 = np.array([ D2[nn*dn:(nn+1)*dn] for nn in self.dofs0 ]) # to [node, D2 flat index 0 to 3]
        return sfcom.mat3d_arr(D2, self.modelplane, reshape=True) # [node,3,3]
        
    def Eij_tranisotropic(self, sn, ei=()):
        """
        Bulk enhancement factors wrt ei=(e1,e2,e3) axes for *transversely isotropic* grains

        If ei=() then CPO eigenframe is used
        """

        nlm = self._nlm_nodal(sn) # [node, component]
        mi, lami = sfcom.eigenframe(nlm, symframe=self.symframe, modelplane=self.modelplane) # [node,i,xyz], [node,i]
        if   len(ei) == 0:          ei = (mi[:,0], mi[:,1], mi[:,2])
        elif len(ei[0].shape) == 1: ei = sfcom.ei_tile(ei, self.numdofs) # ei[i][node,xyz] = (e1[node,xyz], e2, e3)

        if self.enhancementmodel == 'LTS':
                    
            Eij_grain, alpha, n_grain = self.homoparams
            if n_grain != 1: raise ValueError('only n_grain = 1 (linear viscous) is supported')
            if not(0 <= alpha <= 1): raise ValueError('alpha should be between 0 and 1')
            Eij = self.sf.Eij_tranisotropic_arr(nlm, *ei, Eij_grain, alpha, n_grain) # [node, Voigt index 1--6]
            # The enhancement factor model depends on effective (homogenized) grain parameters, calibrated against deformation tests.
            # For CPOs far from the calibration states, negative values *may* occur where Eij should tend to zero if truncation L is not large enough.
            Eij[Eij < 0] = 1e-2 # Set negative E_ij to a very small value (flow inhibiting)
        
        elif self.enhancementmodel == 'APEX':

            Emin, Emax = self.homoparams
            Eij = self.sf.Eij_tranisotropic_APEX_arr(nlm, *ei, Emin, Emax) # [node, Voigt index 1--6]
            
        return self._np2func(ei, Eij, lami)
        
    # @TODO func inp args to use self.homoparams
    def Eij_orthotropic(self, sb, sn, Eij_grain, alpha, n_grain, ei=()):
        """
        Same as Eij_tranisotropic() but for *orthotropic* grains
        """
    
        if n_grain != 1: raise ValueError('only n_grain = 1 (linear viscous) is supported')
        if not(0 <= alpha <= 1): raise ValueError('alpha should be between 0 and 1')
       
        blm = self._nlm_nodal(sb) # [node, component]
        nlm = self._nlm_nodal(sn)
        vlm = 0*nlm # calculate from (blm,nlm) using joint ODF
        mi, lami = sfcom.eigenframe(nlm, symframe=self.symframe, modelplane=self.modelplane) # [node,i,xyz], [node,i]
        if   len(ei) == 0:           ei = (mi[:,0], mi[:,1], mi[:,2])
        elif len(ei[0].shape) == 1:  ei = sfcom.ei_tile(ei, self.numdofs) # ei[i][node,xyz] = (e1[node,xyz], e2, e3)
        Eij = self.sf.Eij_orthotropic_arr(blm, nlm, vlm, *ei, Eij_grain, alpha, n_grain) # [node, Voigt index 1--6]

        Eij[Eij < 0] = 1e-2 # Set negative E_ij to a very small value (flow inhibiting)
        
        return self._np2func(ei, Eij, lami)
        

    def E_CAFFE(self, sn, u, n_grain=1, Emin=0.1, Emax=10):
        """
        CAFFE model (Placidi et al., 2010)
        """
        Df = project( sym(grad(u)), self.G).vector()[:] # flattened strain-rate tensor
        E_CAFFE = Function(self.R)
        E_CAFFE.vector()[:] = self.sf.E_CAFFE_arr(self._nlm_nodal(sn), self.mat3d(Df), Emin, Emax, n_grain)[:]
        return E_CAFFE
        
            
    def E_EIE(self, u, Eij, mi, n, q=None, dim=2):
        """
        Equivalent isotropic enhancement factor (EIE) estimated from tensorial viscous structure.
        """

        if q is None: q = -n-1 # ideal solution
        
        D = sym(grad(u))
        if dim==2: D = as_tensor(sfcom.mat3d(D, self.modelplane)) # 2x2 to 3x3 strain-rate tensor
        
        ort = Orthotropic(n=Constant(n))
        iso = Isotropic(  n=Constant(n))

        # Effective strain rates
        epsE_ort = sqrt(inner(ort.C_inv(D, mi, Eij), D))
        epsE_iso = sqrt(inner(iso.C_inv(D),D))

        # Equivalent isotropic enhancement factor
        E = (epsE_ort/epsE_iso)**(q)
        
        return project(E, self.R) # E as function


    def pfJ(self, s, *args, **kwargs):
        """
        Pole figure J (pfJ) index
        """
        J = Function(self.R)
        J.vector()[:] = sfcom.pfJ(self._nlm_nodal(s), *args, **kwargs)[:]
        return J
        

    def shearfrac_SSA(self, u2):
        """
        Shear fraction for SSA flows (Graham et al., 2018)
        """

        Q = FunctionSpace(self.mesh, 'DG', 0)

        if self.modelplane == 'xy':
            z = as_vector([0,0,1])
            u3 = as_vector([u2[0],u2[1],0])
            
        elif self.modelplane == 'xz':
            z = as_vector([0,1,0])
            u3 = as_vector([u2[0],0,u2[1]])
        
        omghatD = z # if SSA
        q = cross(u3, omghatD) # normal
#        n = project(q/sqrt(dot(q,q)), Qv) # normalize
        n = q/sqrt(dot(q,q)) # normalize
        
        D2 = sym(grad(u2))
        D3 = as_tensor(sfcom.mat3d(D2, self.modelplane)) # 3x3 strain-rate tensor
        
        F = dot(D3,n) - dot(n,dot(D3,n))*n - dot(omghatD,dot(D3,n))*omghatD 
        eps_prime = sqrt(dot(F,F)) # eqn (7) 
        
        eps_E = sqrt(inner(D3,D3)/2)
        gamma = project(eps_prime/eps_E, Q) # shear fraction; normalized by effective strain rate
        
        return (gamma, u3, D3, n, Q)
        
        
    def coaxiality(self, A, B):
        """ 
        Measure of coaxiality between A and B
        """
    
        C = np.matmul(A,B) - np.matmul(B,A) # commutator
        C2 = np.linalg.norm(C) # Frobenius (two) norm: sqrt(C^T : C)
        coaxiality = C2
        return coaxiality
        

    def chi(self, u2, s, verbose=False):
        """
        Fabric compatibility measure \chi
        """

        gam, u3, D3, n, Q = self.shearfrac_SSA(u2)
        G = FunctionSpace(self.mesh, 'DG', 0)
        gam = project(gam, G)

        D_np = project(D3/sqrt(inner(D3,D3)), self.T).vector()[:] # normalized, flattened entries for all nodes
                
        N = outer(n,n) # shear plane normal
        N_np = project(N/sqrt(inner(N,N)), self.T).vector()[:] # normalize, flattened entries for all nodes

        ### Construct chi nodal-wise

        nlm = self._nlm_nodal(s)
        chi_np = np.zeros((self.numdofs))
        gam_np = gam.vector()[:]
        
        for nn in self.dofs0: 
            
            a2 = self.sf.a2(nlm[nn,:])
            a2 /= np.linalg.norm(a2) 

            I3 = np.arange(nn*9,(nn+1)*9) # flattened entries for tensor at node nn
            D = np.reshape(D_np[I3], (3,3))
            N = np.reshape(N_np[I3], (3,3))

#            offset, rescale = 1/2, 2
#            Dc = rescale*((1-self.coaxiality(D, a2)) - offset)
#            
#            offset, rescale = 1-np.sin(np.deg2rad(45)), 1/np.sin(np.deg2rad(45))
#            Nc = rescale*( (1-self.coaxiality(N, a2)) - offset) #np.sin(np.deg2rad(45))
            
            Dc = 1 - 2*self.coaxiality(D, a2)
            Nc = 1 -   self.coaxiality(N, a2)/np.sin(np.deg2rad(45))
            
            chi_np[nn] = np.sqrt((1-gam_np[nn])*Dc**2 + gam_np[nn]*Nc**2)
    
            if verbose:                    
                np.set_printoptions(precision=3)
                print('------------')
                print('---- a2 ----'); print(a2)
                print('---- D ----');  print(D)
                print('---- N ----');  print(N)
                print('chi = %.2f -- coax(a2,D) = %.2f --  coax(a2,N) = %.2f -- gamma = %.2f'%(chi_np[nn], Dc, Nc, gam_np[nn]))

        chi = Function(G)
        chi.vector()[:] = chi_np[:]
        chi = project(chi, Q) # element-based field
        
        return chi
        
