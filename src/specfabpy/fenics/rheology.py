#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022-

#import code # code.interact(local=locals())
import numpy as np
from dolfin import *
from ..common import *


class Isotropic():

    """
    Isotropic viscoplastic rheology
    """

    def __init__(self, n=1):
        self.n = n

    def D(self, S, A, *args, **kwargs): return self.fluidity( S, A, *args, **kwargs) * self.C_fwd(S, A, *args, **kwargs)
    def S(self, D, A, *args, **kwargs): return self.viscosity(D, A, *args, **kwargs) * self.C_inv(D, A, *args, **kwargs)

    def fluidity( self, S, A, *args, **kwargs): return              A * inner(S,S)**((self.n-1)/2) # eta^-1
    def viscosity(self, D, A, *args, **kwargs): return A**(-1/self.n) * inner(D,D)**((1-self.n)/(2*self.n)) # eta
    
    def C_fwd(self, S, *args, **kwargs): return S
    def C_inv(self, D, *args, **kwargs): return D
        
        
class IsotropicPlastic():

    """
    Isotropic plastic rheology
    """

    def __init__(self, *args, **kwargs):
        pass
        
    def D(self, S, *args, **kwargs): return self.fluidity( S, *args, **kwargs) * self.C_fwd(S, *args, **kwargs)
    def S(self, D, *args, **kwargs): return self.viscosity(D, *args, **kwargs) * self.C_inv(D, *args, **kwargs)

    def fluidity( self, S, *args, **kwargs): return None # not supported
    def viscosity(self, D, *args, **kwargs): return inner(D,D)**(-1/2) # for n->inf => second invariant of the strainrate tensor
        
    def C_fwd(self, S, *args, **kwargs): return S
    def C_inv(self, D, *args, **kwargs): return D
        

class Orthotropic():

    """
    Orthotropic viscoplastic rheology
    """

    def __init__(self, n=1):
        self.n = n
        self.ji = [1,2,0]
        self.ki = [2,0,1]

    def D(self, S, A, mi, Eij, *args, **kwargs): return self.fluidity( S, A, mi, Eij, *args, **kwargs) * self.C_fwd(S, mi, Eij, *args, **kwargs)
    def S(self, D, A, mi, Eij, *args, **kwargs): return self.viscosity(D, A, mi, Eij, *args, **kwargs) * self.C_inv(D, mi, Eij, *args, **kwargs)        

    def fluidity(self, S, A, mi, Eij, *args, modelplane=None, **kwargs):
        S3 = S if modelplane is None else as_tensor(mat3d(S, modelplane)) # 2D to 3D with zero components in dimension not resolved
        C3 = self.C_fwd(S3, mi, Eij, *args, modelplane=None, **kwargs)
        etainv = A * inner(C3,S3)**((self.n-1)/2)
        return etainv 
    
    def viscosity(self, D, A, mi, Eij, *args, modelplane=None, **kwargs):
        D3 = D if modelplane is None else as_tensor(mat3d(D, modelplane)) # 2D to 3D with zero components in dimension not resolved
        C3 = self.C_inv(D3, mi, Eij, *args, modelplane=None, **kwargs)
        eta = A**(-1/self.n) * inner(C3,D3)**((1-self.n)/(2*self.n))  
        return eta 
    
    def _mimi(self, mi):
        return [outer(mi[ii],mi[ii]) for ii in range(3)]
       
    def _Mi_shear(self, X, mi):
        Mi = [sym(outer(mi[self.ji[ii]], mi[self.ki[ii]])) for ii in range(3)]
        Ji = [inner(X, Mi[ii]) for ii in range(3)]
        return (Mi, Ji)
    
    def C_fwd(self, S, mi, Eij, *args, modelplane=None, **kwargs):
    
        if modelplane is None: pass
        else: S = as_tensor(mat3d(S, modelplane)) # 2D to 3D with zero components in dimension not resolved
    
        # Tensors and invariants
        mimi = self._mimi(mi)
        Mi = [(mimi[self.ji[ii]] - mimi[self.ki[ii]])/2 for ii in range(3)]
        Ji = [inner(S, Mi[ii]) for ii in range(3)]
        (Mi3, Ji3) = self._Mi_shear(S, mi)

        # Coefficients
        nn = 2/(self.n+1)
        E11, E22, E33, E23, E13, E12 = Eij
        k = 4/3
        etai = [k*(E22**nn + E33**nn - E11**nn), \
                k*(E33**nn + E11**nn - E22**nn), \
                k*(E11**nn + E22**nn - E33**nn)]
        etai3 = [2*E23**nn, 2*E13**nn, 2*E12**nn]
        
        # Construct tensorial part of flow law
        C = sum([ etai[ii]*Ji[ii]*Mi[ii] + etai3[ii]*Ji3[ii]*Mi3[ii] for ii in range(3) ])
        
        if modelplane is None: pass
        else: C = as_tensor(mat2d(C, modelplane)) # 3D back to 2D model plane (disregarding components along dimension not resolved)
        
        return C
        
    def C_inv(self, D, mi, Eij, *args, modelplane=None, **kwargs):
    
        if modelplane is None: pass
        else: D = as_tensor(mat3d(D, modelplane)) # 2D to 3D with zero components in dimension not resolved

        # Tensors and invariants
        mimi = self._mimi(mi)
        Mi = [(Identity(3)-3*mimi[ii])/2 for ii in range(3)]
        Ji = [-3/2*inner(D, mimi[ii]) for ii in range(3)] # = Ji:D since I:D = tr(D) = 0
        (Mi3, Ji3) = self._Mi_shear(D, mi)
  
        # Coefficients
        nn = 2/(self.n+1)
        E11, E22, E33, E23, E13, E12 = Eij
        gamma  = 2*(E22*E33)**nn + 2*(E33*E11)**nn + 2*(E11*E22)**nn - E11**(2*nn) - E22**(2*nn) - E33**(2*nn)         
        k = 4/(3*gamma)
        etai = [k*(E22**nn + E33**nn - E11**nn), \
                k*(E33**nn + E11**nn - E22**nn), \
                k*(E11**nn + E22**nn - E33**nn)]
        etai3 = [2*E23**(-nn), 2*E13**(-nn), 2*E12**(-nn)]
        
        # Construct tensorial part of flow law
        C = sum([ etai[ii]*Ji[ii]*Mi[ii] + etai3[ii]*Ji3[ii]*Mi3[ii] for ii in range(3) ])

        if modelplane is None: pass
        else: C = as_tensor(mat2d(C, modelplane)) # 3D back to 2D model plane (disregarding components along dimension not resolved)

        return C


class ASSA():

    """
    Anisotropic Shallow Shelf Approximation (ASSA)
    """
    
    def __init__(self, n=3, rheology='Orthotropic', modelplane='xz', **kwargs):

        self.n = n
        self.rheology = rheology
        self.modelplane = modelplane
        
        if self.rheology not in ['Isotropic', 'Orthotropic']: 
            raise ValueError('Rheology "%s" not supported. Supported rheologies are "Orthotrpic" or "Isotropic"'%(self.rheology))
        else:
            if self.rheology == 'Isotropic' :   self.flowlaw = Isotropic(n=n)
            if self.rheology == 'Orthotropic' : self.flowlaw = Orthotropic(n=n)
       
    def strainrate(self, u): 
        return sym(grad(u)) 

    def A(self, T, T0=0):
        A0, Q, R = 3.985e-13, 60e3, 8.314 
        A = A0*exp(-Q/(R*(T+T0))) # flow-rate factor 
        return A

    def get_R(self, u, n, A, mi, Eij):
    
        """
        Resistive stress tensor for SSA stress balance
        """
        
        D2 = self.strainrate(u) # 2x2 strain-rate tensor
        D = as_tensor(mat3d(D2, self.modelplane)) # 3x3 strain-rate tensor
        C = self.flowlaw.C_inv(D, mi, Eij)

        eta = A**(-1/n) * (1/2*inner(C,D))**((1-n)/(2*n)) # nonlinear viscosity                                         
        tau = eta * C # deviatoric stress tensor

        i0, j0 = 0,0
        if   self.modelplane=='xy': i1, j1 = 1,1
        elif self.modelplane=='xz': i1, j1 = 2,2
        
        R = as_tensor([ \
            [2*tau[i0,j0] + tau[i1,j1],   tau[i0,j1]             ], \
            [               tau[i1,j0], 2*tau[i1,j1] + tau[i0,j0]] \
        ])

        return R
        
            
