#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024

import numpy as np
from .constants import *
from dolfin import *

class MomentumBalance():

    """
    Momentum balance (Stokes flow BVP)
    """
    
    def __init__(self, mesh, **kwargs):
        
        self.meshobj = mesh
        self.mesh, self.boundaries, self.ds, self.n = mesh.mesh, mesh.boundaries, mesh.ds, mesh.n
        self.ex, self.ey = Constant((1,0)), Constant((0,1))
        self.set_geometry(self.meshobj.h, self.meshobj.s) # set initial thickness 
        self.Eiso = Constant(1)
        
        ### Funtion spaces
        
        self.Uele = VectorElement("CG", self.mesh.ufl_cell(), 2) # velocity element
        self.U  = FunctionSpace(self.mesh, self.Uele)
        self.u  = Function(self.U)      # unknown in NONLINEAR problem
        self.u0 = TrialFunction(self.U) # unknown in LINEAR problem
        self.v  = TestFunction(self.U)  # weight function

        self.Q  = FunctionSpace(self.mesh, "DG", 0) # for A (flow rate factor)

        ### Boundary conditions for stokes problem (assumed constant)
        
        self.bcs = []
        self.bcs += [ DirichletBC(self.U.sub(1), Constant(0), self.boundaries, DOM_ID__BOT) ];
        self.bcs += [ DirichletBC(self.U.sub(1), Constant(0), self.boundaries, DOM_ID__TOP) ];
        self.bcs += [ DirichletBC(self.U.sub(0), Constant(0), self.boundaries, DOM_ID__TOP) ];
                            
    def set_rheology(self, rheology): 
        self.rheology = rheology 
        if self.rheology.n==2: Eadj = 1e5
        if self.rheology.n==3: Eadj = 1
        if self.rheology.n==4: Eadj = 8e-5
        self.A0 = interpolate(Constant(Eadj * self.rheology.A(self.rheology.T)), self.Q) # done this way only because temperature is assumed constant 
        
    def set_fabric(self, fabric): self.fabric = fabric
    def set_geometry(self, h, s): self.h, self.s = h, s

    def set_E_EIE(self, *args):
        self.Eiso = self.fabric.enhancementfactor.E_EIE(*args)
        return self.Eiso

    def set_E_CAFFE(self, *args):
        self.Eiso = self.fabric.enhancementfactor.E_CAFFE(*args)
        return self.Eiso
            
    def solve(self, u_guess, tolmul=1, relaxation=0.45, maxiter=150):

        # Linear solve?
        if (self.rheology.n == 1) or (u_guess is None):
            A = self.Eiso * self.A0 
            F = self.get_functional(self.u0,self.v, A, 1)
            u_lin = Function(self.U)
            solve(lhs(F)==rhs(F), u_lin, self.bcs) 
            self.u.assign(u_lin) # linear solution
        else:
            u_guess.set_allow_extrapolation(True) 
            LagrangeInterpolator.interpolate(self.u, u_guess) # init guess for nonlin solver was passed directly

        # Nonlinear solve? Then use the init guess stored in "u"
        if self.rheology.n > 1:   
            A = self.Eiso * self.A0 # Apply isotropic enhancement factor  
            F = self.get_functional(self.u,self.v, A, self.rheology.n) 
            kwargs = {"relative_tolerance": tolmul*1e-4, "absolute_tolerance": tolmul*1e-3, "relaxation_parameter":relaxation, \
                       "maximum_iterations": maxiter, 'linear_solver': 'gmres', 'preconditioner': 'ilu'}
            solve(F == 0, self.u, self.bcs, solver_parameters={"newton_solver": kwargs})
            

    def get_functional(self, u,v, A, n):

        """
        Stokes functional
        """

        # Driving stress        
        F = Constant(DENSITY_ICE*GRAV_ACCEL)*self.h*inner(grad(self.s),v) * dx 
        
        # Viscous stresses
        R = self.rheology.get_R(u, n, A, self.fabric.mi, self.fabric.Eij)
        F += inner(self.h*R, sym(grad(v))) * dx

        # Pressure/stress BCs
        F += -self.meshobj.p_out*inner(self.n, v)*self.ds(DOM_ID__OUT)

        return F

