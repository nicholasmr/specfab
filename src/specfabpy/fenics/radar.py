#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
FEniCS interface for GTM radar model (Rathmann et al., 2022)
"""

import copy, code # code.interact(local=locals())

import numpy as np

from .ice import *
from ..tamm.layerstack import LayerStack

from scipy.interpolate import LinearNDInterpolator
import dolfin as df


class IceRadar(IceFabric):
    
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
    

    def build_radar_transect(self, 
            x, zb, zs, Nz=50, 
            freq=300e6, epsa=3.16, deps=0.034, mu=1, modeltype='GTM'
        ):

        ### Geometry 

        self.nx = len(x) # Number of ice columns to build radar matrix model on    
        self.nz = Nz     # Number of vertical layers in radar model (irrespective of ice thickness)        
                
        self.x  = x  # Location of radar columns
        self.zs = zs # Surface interface
        self.zb = zb # Basal interface

        ### Wave propagation constants 

        self.freq = freq # frequency [Hz]
        self.omega = 2*np.pi*self.freq # angular frequency

        ### Initialize transfer matrix model per column (per x location)
        
        # The generalized matrix model (Rathmann et al., 2022) is used to also support fabric fields *without* a vertical eigenvector

        kw_stack = dict(epsa=[epsa], epsc=[epsa+deps], f=self.freq, modeltype=modeltype)

        if not hasattr(self, 'layerstack'): # stack already set?
            
            self.X = np.zeros((self.nx,self.nz))
            self.Z = np.zeros((self.nx,self.nz))
            
            self.stack_para = []
            self.stack_perp = []
            
            zipped = zip(self.x, self.zs, self.zb)
            
            for ii,(xi,zsi,zbi) in enumerate(zipped):
            
                zj = np.linspace(zsi, zbi, self.nz)
                self.X[ii,:], self.Z[ii,:] = xi, zj
                
                nlm = np.zeros((len(zj),self.nlm_len_full), dtype=np.complex64)
                
                for jj, zjj in enumerate(zj):
                
                    nlm[jj,:] = self.get_nlm(xi,zjj)
                    # surface layer must be isotropic; fix small numerical errors on order of 1e-20
                    if jj==0: 
                        nlm[jj,1:] = 0 
                
                # Build two matrix models per radar column: one for HH and one for VV acquisitions 
                # corresponding to horizontal in- and out-of-model plane directions
                ls_para = LayerStack(nlm, zj, beta=[0,],              **kw_stack)
                ls_perp = LayerStack(nlm, zj, beta=[np.deg2rad(90),], **kw_stack)
                
                self.stack_para.append(ls_para)
                self.stack_perp.append(ls_perp)
                
       
    def np2fs(self, X, Z, V):
    
        ### Make functions (CG1) based on data in numpy array

        Qele = ('CG', 1) # for projecting scalar fields
        Q = FunctionSpace(self.mesh, *Qele)

        points = np.column_stack([X.ravel(), Z.ravel()])
        f_interp = LinearNDInterpolator(points, V.ravel()) # V = vals
        f = df.Function(Q)
        gdim = self.mesh.geometry().dim()
        coords = Q.tabulate_dof_coordinates().reshape((-1, gdim))
        f.vector()[:] = f_interp(coords)
        
        return f
        
        
    def get_Si(self):

        ### Slownesses of nadir propagating waves

        S1_np,S2_np = np.zeros((self.nx,self.nz)), np.zeros((self.nx,self.nz))
        for ii in range(self.nx):
            ki = np.real(self.stack_para[ii].get_propconst())[:,[0,1]]
            S1_np[ii,:] = ki[:,0]/self.omega
            S2_np[ii,:] = ki[:,1]/self.omega
            
        S1 = self.np2fs(self.X, self.Z, S1_np)
        S2 = self.np2fs(self.X, self.Z, S2_np)

        return (S1,S2)
        

    def get_dTWTT(self, ascycles=False, BOUNDARY_TOP=1):

        ### Two-Way Travel Time (TWTT) difference of horizontally eigenpolarized waves traveling to some depth z, and back

        V = df.FunctionSpace(self.mesh, "CG", 1)
        w = df.TestFunction(V)
        dTWTT = df.TrialFunction(V)
        
        S1, S2 = self.get_Si() # Use transfer matrix model to determine (sorted) wave velocities
        dSi = S1-S2 
        a = df.Dx(dTWTT,1)*w*df.dx
        L = 2*dSi*w*df.dx

        dTWTT_sol = df.Function(V)
        bc = df.DirichletBC(V, df.Constant(0.0), self.boundaries, BOUNDARY_TOP)
        df.solve(a == L, dTWTT_sol, bc)

        dphi = df.Function(V)
        vals = self.omega * dTWTT_sol.vector()[:]
#        if wrap:     vals = vals % 2*np.pi # wrapped phase into [0, 2π)
        if ascycles: vals /= 2*np.pi       # no. of cycles instead
        
        dphi.vector()[:] = vals

        return dTWTT_sol, dphi


    def get_Pij(self, E0=1):

        ### HH and VV return power 

        Phh_np, Pvv_np = np.zeros((self.nx,self.nz)), np.zeros((self.nx,self.nz))
        
        for ii in range(self.nx):
            Phh_np[ii,:-1], *_ = self.stack_para[ii].get_returns(E0=E0)
            Pvv_np[ii,:-1], *_ = self.stack_perp[ii].get_returns(E0=E0)
            Phh_np[ii,-1], Pvv_np[ii,-1] = Phh_np[ii,-2], Pvv_np[ii,-2]

        Phh = self.np2fs(self.X, self.Z, Phh_np)
        Pvv = self.np2fs(self.X, self.Z, Pvv_np)
        
        return (Phh,Pvv)


    def get_n2m(self):
    
        ### Return absolute value of n20, n21, n22 coefficients on DG0
    
        nlm_np = self.enhancementfactor._nlm_nodal(self.s) 

        n20 = Function(self.enhancementfactor.R)
        n21 = Function(self.enhancementfactor.R)
        n22 = Function(self.enhancementfactor.R)

        n20.vector()[:] = np.abs(nlm_np[:,self.sf.I20+0]) 
        n21.vector()[:] = np.abs(nlm_np[:,self.sf.I20+1]) 
        n22.vector()[:] = np.abs(nlm_np[:,self.sf.I20+2]) 

        return n20, n21, n22
        
    
#    def get_dlam(self):
#    
#        ### Delta lambda and lambda_z 
#    
#        nlm = self.enhancementfactor._nlm_nodal(self.s) 
#        mi, lami = sfcom.eigenframe(nlm, modelplane='dl')
#        mi_df, _, lami_df = self.enhancementfactor._np2func( (mi[:,0], mi[:,1], mi[:,2]), None, lami )
#        dlam_df = lami_df[0]-lami_df[1] # horizontal eigenvalue difference, where dlam>0 means that the x eigenvalue is larger than the y eigenvalue
#        lamz_df = lami_df[2]
#        return dlam_df, lamz_df

