#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024

import numpy as np
from datetime import datetime
import os, pickle

from .constants import *
#from .defaultsettings import * 

from .mesh import *
from .fabric import *
from .rheology import *
from .momentum import *

from dolfin import *

class FlowModel():

    """
    Flow model
    """
    
    def __init__(self, Nt, verbose=True, kwargs_mesh={}, kwargs_rheo={}, kwargs_fabric={}):
    
        """ 
        Initialize model components
        
        kwargs_* values are used to override default values in defaultsettings.py 
        
        Inputs:
        Nt      Total number of integration steps to be taken
        """

        self.Nt = Nt 
        self.tt   = 0 # total number of integration steps
        self.time = 0 # total integration time (seconds)
        self.timeax = np.full((self.Nt+1),np.nan) # [0, dt_0, dt_0+dt_1, ...]
        self.timeax[0] = 0
        
        ### MPI
        self.comm = MPI.comm_world
        self.rank = MPI.rank(self.comm)
        
        ### Mesh
        self.mesh = MeshGeom(**kwargs_mesh)
        
        ### Fabric and rheology
        fabric        = Fabric(self.mesh, **kwargs_fabric)
        self.rheology = Rheology(**kwargs_rheo)

        ### Momentum balance (Stokes BVP)
        self.momentum = MomentumBalance(self.mesh)
        self.momentum.set_rheology(self.rheology) 
        self.momentum.set_fabric(fabric) 
        
        ### Print setup
        if verbose:
            from tabulate import tabulate
            print(tabulate([['Rheology', self.rheology.rheology], ['nglen', self.rheology.n], ['T (deg. C)', self.rheology.T-273.15], \
                            ['Fabric dynamics', self.momentum.fabric.fabdyn], ['Fabric symm. frame', self.momentum.fabric.symframe], \
                            ['Fabric model plane', self.momentum.fabric.modelplane], ['L', self.momentum.fabric.L], \
                            ['alpha', self.momentum.fabric.alpha], ['Eij_grain', self.momentum.fabric.Eij_grain]
                            ], headers=['Parameter', 'Value'], tablefmt="rst"))
            print()

    def integrate(self, dt, u_guess=None, tolmul=1e-1, relaxation=0.25, maxiter=250, verbose=True):
        
        """
        Integrates model by dt (time step length)
        """
        
        tstep = datetime.now() # for timing how long integration took
        
        ### Solve Stokes problem
        MPI.barrier(self.comm)                
        if verbose: self.print_statusmsg_start('Solving momentum balance (tolmul=%.2e, relaxation=%.2f, maxiter=%i)'%(tolmul, relaxation, maxiter))
        self.momentum.solve(u_guess if self.tt==0 else self.momentum.u, tolmul=tolmul, relaxation=relaxation, maxiter=maxiter) # updates momentum.u (velocity field)
        self.print_statusmsg_end()
       
        ### Solve fabric evolution
        if self.rheology.rheology != 'Isotropic' and dt>0:
            MPI.barrier(self.comm)
            self.print_statusmsg_start('Solving fabric evolution')
            S = sym(grad(self.momentum.u)) # assume strain-rate tensor is approximately co-axial to stress tensor
            self.momentum.fabric.evolve(self.momentum.u, S, self.rheology.T, dt) # update fabric state vector
            self.print_statusmsg_end()
            
        ### Solve thickness evolution
        if dt>0:
            self.mesh.update(self.momentum.u, dt)
            self.momentum.set_geometry(self.mesh.h, self.mesh.s)
            
        ### Update counters
        self.tt   += 1
        self.time += dt
        self.timeax[self.tt] = self.time
        
        if verbose: 
            print("[OK] Integration step completed, took %.1f seconds\n"%(datetime.now()-tstep).total_seconds())
        
        
    def dt_CFL(self, C_CFL=0.5):
        h_min = self.mesh.mesh.hmin()
        v_max = abs(self.momentum.u.vector()[:]).max()
        dt = C_CFL*h_min/v_max
        return dt
        
        
    def print_statusmsg_start(self, string):
        print('[->] %s'%(string))
        self.tstart=datetime.now()
        
                
    def print_statusmsg_end(self):
        if self.tstart is not None: 
            print("[OK] Took %.1f seconds\n"%(datetime.now()-self.tstart).total_seconds())
        self.tstart = None
        
   
    def save(self, fname):
        hdf5 = HDF5File(self.mesh.mesh.mpi_comm(), fname, "w")
        hdf5.write(self.mesh.mesh, "mesh")
        hdf5.write(self.momentum.h, "h")
        hdf5.write(self.momentum.u, "u")
        hdf5.write(self.momentum.A0, "A0")
        hdf5.write(self.momentum.fabric.s, "s") # fabric state vector        
        hdf5.close()

    def set_state(self, u, h, A0, s):
        self.mesh.set_thickness(h)
        self.momentum.u.assign(u)
        self.momentum.set_geometry(self.mesh.h, self.mesh.s)
        self.momentum.A0 = A0
        self.momentum.fabric.set_state(s) # runs update_Eij() automatically


def load_state(fname, L, mesh=None):

    print('[++] Reloading saved state %s'%(fname))
    
    if mesh is None:
        mesh = Mesh()
        hdf5 = HDF5File(mesh.mpi_comm(), fname, "r")
        hdf5.read(mesh, '/mesh', False)
    else:
        hdf5 = HDF5File(mesh.mpi_comm(), fname, "r")
    
    H = FunctionSpace(mesh, "CG", 1)
    h = Function(H)
    hdf5.read(h, "/h")
    
    U  = FunctionSpace(mesh, VectorElement("CG", mesh.ufl_cell(), 2))
    u  = Function(U)        
    hdf5.read(u, "/u")

    Q = FunctionSpace(mesh, "DG", 0)
    A0 = Function(Q)
    hdf5.read(A0, "/A0")

    cpo = CPO(mesh, None, L)
    s = Function(cpo.V)
    hdf5.read(s, "/s")

    return (u, h, A0, s)
    
