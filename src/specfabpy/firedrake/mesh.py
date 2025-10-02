#!/usr/bin/python3
# Nicholas Rathmann, 2025-

import firedrake as fd
import numpy as np
from scipy import interpolate
#import copy, sys, time, code # code.interact(local=locals())

class PeriodicBumpyMesh():

    """
    xz mesh with cosine bumpy bed and periodic boundaries
    """

    name = 'PeriodicBumpyMesh'
    nameshort = 'bumpybed'

    def __init__(self, L=20e3, f_H0=None, Abump=250, Nbump=1, angle=10, nx=50, nz=25, fsreg=1e-6):
        
        ### Save inputs
        
        self.L     = L     # domain length
        self.Abump = Abump # bump aplitude
        self.Nbump = Nbump # number of bumps
        self.nx    = nx    # horizontal resolution
        self.nz    = nz    # vertical resolution
        self.angle = angle # mesh angle of inclination
        self.fsreg = fsreg # regularization magnitude of free surface problem

        # AUX
        self.Ndofs  = self.nx-1
        self.beta   = np.deg2rad(angle)
        self.hx_min = self.L/self.nx
                
        ### Setup mesh

        self.xs = np.linspace(0, self.L, self.Ndofs) # x-coords of nodes and CG1 DOFs
        zs = f_H0(self.xs/self.L) if f_H0 is not None else 3e3*np.ones((self.Ndofs)) # initial surface profile
        self.update_mesh(self.xs, zs)
        
        self.boundaries = (1,2)
        self.id_bot, self.id_top = self.boundaries # unpack
                
        ### Function spaces for free surface problem
        
        self.smesh = fd.PeriodicIntervalMesh(self.Ndofs, self.L) # mesh for surface height problem
        self.S  = fd.FunctionSpace(self.smesh, "CG", 1)
        self.v  = fd.TestFunction(self.S)  # weight funcs
        self.s  = fd.TrialFunction(self.S) # solution
        self.s0 = fd.Function(self.S)      # previous solution
        self.s0.dat.data[:] = zs           # initial surface profile
               

    def print_setup(self):
        from tabulate import tabulate
        from termcolor import colored
        print()
        cbold = lambda s: colored(s, attrs=['bold'])
        cnum = lambda s: colored(s,'magenta', attrs=["dark"])
        cflg = lambda s: colored(s,'green', attrs=["dark"])
        print(tabulate([['name', self.name], ['L (km)', '%.1f'%(self.L/1e3)], ['Abump (m)', self.Abump], ['Nbump', self.Nbump], \
                        ['nx, nz', (self.nx, self.nz)], ['angle (deg)', self.angle], \
                        [cnum('numerics: fsreg'), "%.1e"%self.fsreg],  \
                        ], headers=[cbold('Mesh parameter'), cbold('Value')], tablefmt="rst"))


    def dt_CFL(self, u):
        
        ux_max = abs(u.sub(0).vector()[:]).max()
        dt_CFL = 0.5*self.hx_min/ux_max
        return dt_CFL # CFL timestep size
                
                
    def update_mesh(self, x_1d, zs_1d):

        self.mesh = fd.PeriodicRectangleMesh(self.nx, self.nz, 1, 1, direction="x")
        x = self.mesh.coordinates.dat.data[:, 0]
        z = self.mesh.coordinates.dat.data[:, 1]
        x0, x1 = min(x), max(x)

        ### Stretch mesh according to basal and surface height profile

        # Bed
        x_new = x0 + x*self.L
        f_zb = lambda x_: self.Abump*(1-np.cos(2*np.pi*self.Nbump*x_))
        zb = f_zb(x)
        
        # Surface 
        f_zs = interpolate.interp1d(x_1d, zs_1d)
        zs = f_zs(x_new)
        z = np.power(z, 1.5) # higher resolution near bed if exponent > 1
        z_new = zb + z*(zs - zb)
        
        # Displace mesh coordinates
        xz_new = np.array([x_new, z_new]).transpose()
        self.mesh.coordinates.dat.data[:] = xz_new 

        # save 1d profiles 
        self.zs = zs_1d
        self.zb = f_zb(x_1d/self.L)
           
           
    def update_surface(self, u, dt): # eps=1e-6 sufficient?
        
        ### Solves free surface equation
        
        usx, usz = self.get_us(u)
        dt_ = fd.Constant(dt)
        a  = self.s *self.v*fd.dx + dt_*usx*self.s.dx(0)*self.v*fd.dx # recall .dx(0) = d/d_x
        L  = self.s0*self.v*fd.dx + dt_*usz*self.v*fd.dx
        a += dt_*fd.Constant(self.fsreg)*fd.inner(fd.grad(self.s), fd.grad(self.v))*fd.dx # regularization
        s_sol = fd.Function(self.S)
        fd.solve(a == L, s_sol, [])
        self.s0.assign(s_sol)
     
        ### Update mesh with new free surface
     
        self.mesh0 = self.mesh # copy of old mesh
        self.update_mesh(self.xs, self.s0.dat.data[:])
     

    def get_us(self, u):
    
        ### Get surface velocities along top (surface) boundary
    
        usx_arr, usz_arr = np.zeros(self.Ndofs), np.zeros(self.Ndofs)
        for ii, x in enumerate(self.xs):
            z = self.s0.dat.data[ii]
            usx_arr[ii], usz_arr[ii] = u((x,z))

        usx, usz = fd.Function(self.S), fd.Function(self.S)
        usx.dat.data[:] = usx_arr[:]
        usz.dat.data[:] = usz_arr[:]
                
        return usx, usz


class HalfDivideMesh():

    """
    xz mesh of a half-width ice divide
    """

    name = 'HalfDivideMesh'
    nameshort = 'halfdivide'

    def __init__(self, L=20e3, f_H0=None, Abump=100, Nbump=1, nx=50, nz=25, fsreg=1e-6):
        
        ### Save inputs
        
        self.L     = L     # domain length
        self.nx    = nx    # horizontal resolution
        self.nz    = nz    # vertical resolution
        self.Abump = Abump # bump aplitude
        self.Nbump = Nbump # number of bumps
        self.fsreg = fsreg # regularization magnitude of free surface problem

        self.beta=self.angle = 0 # mesh slope
        
        # AUX
        self.Ndofs  = self.nx-1
        self.hx_min = self.L/self.nx
                
        ### Setup mesh

        self.xs = np.linspace(0, self.L, self.Ndofs) # x-coords of nodes and CG1 DOFs
        zs = f_H0(self.xs/self.L) if f_H0 is not None else 3e3*np.ones((self.Ndofs)) # initial surface profile
        self.update_mesh(self.xs, zs)
        
        self.boundaries = (1,2,3,4)
        self.id_lft, self.id_rht, self.id_bot, self.id_top = self.boundaries # unpack
                
#        ### Function spaces for free surface problem
#        
#        self.smesh = fd.IntervalMesh(self.Ndofs, self.L) # mesh for surface height problem
#        self.S  = fd.FunctionSpace(self.smesh, "CG", 1)
#        self.v  = fd.TestFunction(self.S)  # weight funcs
#        self.s  = fd.TrialFunction(self.S) # solution
#        self.s0 = fd.Function(self.S)      # previous solution
#        self.s0.dat.data[:] = zs           # initial surface profile
               

    def print_setup(self):
        from tabulate import tabulate
        from termcolor import colored
        print()
        cbold = lambda s: colored(s, attrs=['bold'])
        cnum = lambda s: colored(s,'magenta', attrs=["dark"])
        cflg = lambda s: colored(s,'green', attrs=["dark"])
        print(tabulate([['name', self.name], ['L (km)', '%.1f'%(self.L/1e3)], \
                        ['nx, nz', (self.nx, self.nz)], \
                        [cnum('numerics: fsreg'), "%.1e"%self.fsreg],  \
                        ], headers=[cbold('Mesh parameter'), cbold('Value')], tablefmt="rst"))


    def dt_CFL(self, u):
        
        ux_max = abs(u.sub(0).vector()[:]).max()
        dt_CFL = 0.5*self.hx_min/ux_max
        return dt_CFL # CFL timestep size
                
                
    def update_mesh(self, x_1d, zs_1d):

        self.mesh = fd.RectangleMesh(self.nx, self.nz, 1, 1)
        x = self.mesh.coordinates.dat.data[:, 0]
        z = self.mesh.coordinates.dat.data[:, 1]
        x0, x1 = min(x), max(x)

        ### Stretch mesh according to basal and surface height profile

        # Bed
        x_new = x0 + x*self.L
        f_zb = lambda x_: self.Abump*(1-np.cos(2*np.pi*self.Nbump*x_))
        zb = f_zb(x)
        
        # Surface 
        f_zs = interpolate.interp1d(x_1d, zs_1d)
        zs = f_zs(x_new)
        z = np.power(z, 1.35) # higher resolution near bed if exponent > 1
        z_new = zb + z*(zs - zb)
        
        # Displace mesh coordinates
        xz_new = np.array([x_new, z_new]).transpose()
        self.mesh.coordinates.dat.data[:] = xz_new 

        # save 1d profiles 
        self.zs = zs_1d
        self.zb = f_zb(x_1d/self.L)
           
           
    def update_surface(self, u, dt): # eps=1e-6 sufficient?
        print('>>> update_surface() not implemented for mesh "HalfDivide"')
        pass
        
#        ### Solves free surface equation
#        
#        usx, usz = self.get_us(u)
#        dt_ = fd.Constant(dt)
#        a  = self.s *self.v*fd.dx + dt_*usx*self.s.dx(0)*self.v*fd.dx # recall .dx(0) = d/d_x
#        L  = self.s0*self.v*fd.dx + dt_*usz*self.v*fd.dx
#        a += dt_*fd.Constant(self.fsreg)*fd.inner(fd.grad(self.s), fd.grad(self.v))*fd.dx # regularization
#        s_sol = fd.Function(self.S)
#        fd.solve(a == L, s_sol, [])
#        self.s0.assign(s_sol)
#     
#        ### Update mesh with new free surface
#     
#        self.mesh0 = self.mesh # copy of old mesh
#        self.update_mesh(self.xs, self.s0.dat.data[:])
#     

#    def get_us(self, u):
#    
#        ### Get surface velocities along top (surface) boundary
#    
#        usx_arr, usz_arr = np.zeros(self.Ndofs), np.zeros(self.Ndofs)
#        for ii, x in enumerate(self.xs):
#            z = self.s0.dat.data[ii]
#            usx_arr[ii], usz_arr[ii] = u((x,z))

#        usx, usz = fd.Function(self.S), fd.Function(self.S)
#        usx.dat.data[:] = usx_arr[:]
#        usz.dat.data[:] = usz_arr[:]
#                
#        return usx, usz

