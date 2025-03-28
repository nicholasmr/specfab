#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024

import numpy as np
import copy 
from .constants import *
import matplotlib.tri as tri
from dolfin import *

class MeshGeom():

    """
    Mesh
    """
    
    h_prev = None    
    
    def __init__(self, mesh_overload=None, *args, **kwargs):

        ### Save arguments
        self.__dict__.update(kwargs)

        ### FEM mesh and function space
        if mesh_overload is None:
            try:    self.mesh = Mesh("mesh/%s.xml"%(self.geom)) 
            except: raise ValueError('Invalid mesh geometry "%s"'%(self.geom))
        else:
            self.mesh = mesh_overload
#            self.boundaries = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1, 0)

        self.boundaries = MeshFunction("size_t", self.mesh, "mesh/%s_facet_region.xml"%(self.geom))
        self.n = FacetNormal(self.mesh) # used as boundary normal vector
        self.ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
                
        self.V = FunctionSpace(self.mesh, "CG", 1)
        self.Q = FunctionSpace(self.mesh, "CG", 1) # for vertex_vals_on_Q() only

        ### Coords
        self.coords, self.triang, self.meshpts = self.get_coords()
        self.L = np.amax(self.coords[0,:])
        self.W = np.amax(self.coords[1,:])
        
        ### Initial geometry
        if self.h0profile == 'uniform':    self.set_thickness(Constant(self.H0)) # flat initial ice shelf
        if self.h0profile == 'triangular': self.set_thickness(Expression('H0 + (H1-H0)*(x[0])/L', H0=self.H0, H1=250, L=self.L, degree=1)) # flat initial ice shelf
       
    def get_coords(self, scale=1):
        coords = scale*copy.deepcopy(self.mesh.coordinates().reshape((-1, 2)).T)
        triang = tri.Triangulation(*coords, triangles=self.mesh.cells()) # calculate in this outer scope since it is a shared prerequsite
        meshpts = (coords[0,:], coords[1,:])
        return (coords, triang, meshpts)
        
    def vertexvals_on_Q(self, F_df):
        return project(F_df, self.Q).compute_vertex_values(self.mesh)
        
    def set_thickness(self, h):
        self.h = project(h, self.V) 
        self.update_aux()
        
    def update_aux(self):
        freeboardpct = (1-DENSITY_ICE/DENSITY_OCEAN) # fraction of ice thickness above sea level
        self.S0 = freeboardpct*self.H0
        self.s = project(Constant(freeboardpct)*self.h, self.V)
        h_submerged = (1-freeboardpct)*self.h # ice thickness below sea level at ocean front
        self.p_out =   Constant(0.5*DENSITY_ICE*GRAV_ACCEL)*self.h**2 \
                     - Constant(0.5*DENSITY_OCEAN*GRAV_ACCEL)*h_submerged**2
        self.p_in  = -Constant(0.5*DENSITY_ICE*GRAV_ACCEL)*self.h**2 # constant
        
    def update(self, u, dt, meltrate0=0.00, regmag=2e-4): 
        
        """
        Solve for ice-thickness evolution: dh/dt + div(h*u) = meltrate
        """
        
        if not self.evolve: return # don't updated if user requests constant thickness profile
        
        self.h_prev = self.h.copy()
        
        ### Function space
        v,h = TestFunction(self.V), TrialFunction(self.V) # weight function, unknown
        hprev = Function(self.V) # prev solution
        hprev.assign(self.h)
        
        ### Problem weak form
        dt = Constant(dt)
        a  = h*v*dx + dt*div(u*h)*v*dx 
        a += dt*Constant(regmag)*inner(grad(h), grad(v))*dx # regularization
        meltrate = self.h_prev/self.H0 * Constant(-1 * MYR_TO_MS * meltrate0) # m ice eq. / s (input is in units of per year)
        L  = hprev*v*dx + dt*meltrate*v*dx 

        ### BCs
        bcs = [DirichletBC(self.V, Constant(self.H0), self.boundaries, DOM_ID__IN), ]

        ### Solve
        hnew = Function(self.V) # solution container
        solve(a == L, hnew, bcs)
        self.h.assign(hnew) # save solution
        self.update_aux()
        
