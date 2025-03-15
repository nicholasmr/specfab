#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

import numpy as np
from dolfin import *
from .. import common as sfcom

class FSE():
    
    """
    Finite Strain Ellipse (FSE) solver
    """

    def __init__(self, mesh, boundaries, F0=np.eye(2), ds=None, nvec=None): 

        ### Save inputs
        self.mesh, self.boundaries = mesh, boundaries
        self.ds = ds   if ds   is not None else Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        self.n  = nvec if nvec is not None else FacetNormal(self.mesh)
#        self.mesh, self.boundaries, self.ds, self.n = mesh, boundaries, ds, n
        #self.F0 = F0 # not used

        ### Function spaces
        eletype, eleorder = 'CG', 1
        self.V = TensorFunctionSpace(self.mesh, eletype, eleorder)
        self.p = TrialFunction(self.V)
        self.q = TestFunction(self.V) # weight functions

        ### Solution containers
        self.F      = Function(self.V) # Current solution
        self.F_prev = Function(self.V) # Previous solution
        
        ### Initialize
        I = Constant( ((1,0),(0,1)) ) # init FSE field as isotropic 
        self.F = interpolate(I, self.V)

        ### Boundary conditions for fabric problem (assumed time constant)
        self.bcs = []


    def evolve(self, usol, dt):
        self.F_prev.assign(self.F) # Notice current state must be set!
        L = grad(usol)
        dtinv = Constant(1/dt)    
        f = Constant(0.5) # f=[0;1], f=0.5 => central differencing
        a = inner(dtinv*self.p - (1-f)*dot(L,self.p), self.q) * dx 
        a += inner(dot(usol, nabla_grad(self.p)), self.q) * dx # advection
        L = inner(dtinv*self.F_prev + f*dot(L,self.F_prev), self.q) * dx
        solve(a==L, self.F, self.bcs)
        
    def eigenframe(self, x, y, modelplane=None):
        F = np.reshape(self.F(x,y), (2,2))
#        return sfcom.eigenframe(sfcom.F2C(F), modelplane=modelplane)
        lami, ei = np.linalg.eig(sfcom.F2C(F))
        I = np.flip(lami.argsort()) # largest eigenvalue pair is first entry
        ei, lami = ei[:,I], lami[I]
        return ei.T, lami
#        sr = 1/np.sqrt(eigvals) # stretching ratios

