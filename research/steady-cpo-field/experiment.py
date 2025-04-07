#!/usr/bin/python3
# Nicholas Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Experiment definition class for steady SSA fabric solver
"""

import os, sys, copy, code # code.interact(local=locals())
import numpy as np
import xarray, pickle, pyproj # pyproj needed!
from dolfin import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.tri as tri

import specfabpy.fenics.tools as sffenicstools
from specfabpy.fenics.ice import IceFabric

from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=11)

from tabulate import tabulate


class Experiment():

    ### Input data files

    fvel_default = "~/ice-velocity-maps/antarctica_ice_velocity_450m_v2.nc" # measures velocities
    fbed_default = "~/ice-velocity-maps/BedMachineAntarctica-v3.nc" # bedmachine

    uxname, uyname = 'VX', 'VY'
    xuname, yuname = 'x', 'y'

    ### specfab setup

    L             = 8      # spectral truncation
    nu_multiplier = 0.8    # adjust S^2 regularization by this multiplicative factor (should not be set lower than ~0.8)
    nu_realspace  = 0.2e-3 # free regularization parameter, controls how large spatial gradients in CPO field are permitted

    modelplane      = 'xz' # xz is much faster than xy (renders state vector real-valued)
    BOUNDARY_INFLOW = 1    # subdomain id

    CAFFE_params = (0.1, 10) # Emin, Emax of CAFFE

    ### Mesh
    
    exname  = None   # experiment name
    meshgeo = 'mesh' # .geo mesh file name

    ### Other
    
    velscale = 3.17098e+8 # m/s to m/yr
    lenscale = 1e-3 # m to km
    
    sx = 9 # default diagnostic figure size 
    sy = 7

    def __init__(self, exname, p0=(0,0), p1=(0,0), fvel=None, fbed=None, meshname='mesh'):
    
        self.exname = exname    
        
        self.x0, self.y0 = p0
        self.x1, self.y1 = p1
        self.xlims = [self.x0,self.x1]
        self.ylims = [self.y0,self.y1]
        
        self.x0km, self.y0km = 1e-3*self.x0, 1e-3*self.y0
        self.x1km, self.y1km = 1e-3*self.x1, 1e-3*self.y1
        
        self.fvel = self.fvel_default if fvel is None else fvel
        self.fbed = self.fbed_default if fbed is None else fbed
        self.meshname = meshname

    def preprocess(self, dxy_h=1, dxy_u=1):
    
        """
        Slice and subsample maps of velocity and geometry.
        Results are saved to .h5 and .pkl files.
        """
    
        self.dxy_h = dxy_h
        self.dxy_u = dxy_u
    
        print('*** Preprocessing requested')
        print(tabulate([['exname', self.exname], ['dxy_u', self.dxy_u], ['dxy_h', self.dxy_h], \
                        ['fvel', self.fvel], ['fbed', self.fbed]], \
                        headers=['Parameter', 'Value'], tablefmt="rst"))

        ### Make mesh

        print('*** Generating FE mesh')
        os.system('gmsh -2 -format msh2 %s/%s.geo'%(self.exname,self.meshname))
        os.system('dolfin-convert %s/%s.msh %s/%s.xml'%(self.exname,self.meshname,self.exname,self.meshname))

        ### Load mesh
        
        mesh,boundaries, Q,Q2,V, coords,cells = self.get_mesh()

        # xarray objects for interpolation points (used below)
        coordsQ = Q.tabulate_dof_coordinates().reshape((-1, 2)).T
        xQ = xarray.DataArray(coordsQ[0,:], dims="z")
        yQ = xarray.DataArray(coordsQ[1,:], dims="z")
        
        coordsQ2 = Q2.tabulate_dof_coordinates().reshape((-1, 2)).T
        xQ2 = xarray.DataArray(coordsQ2[0,:], dims="z")
        yQ2 = xarray.DataArray(coordsQ2[1,:], dims="z")

        ### Slice and coarsen

        print('*** Subsampling .nc fields')

        ds_vel0 = xarray.open_mfdataset(self.fvel)
        ds_bed0 = xarray.open_mfdataset(self.fbed)
        xrng = slice(self.x0, self.x1)
        yrng = slice(self.y1, self.y0)
        ds_u = ds_vel0.sel(x=xrng,y=yrng).coarsen(x=self.dxy_u, boundary='trim').mean().coarsen(y=self.dxy_u, boundary='trim').mean()
        ds_b = ds_bed0.sel(x=xrng,y=yrng).coarsen(x=self.dxy_h, boundary='trim').mean().coarsen(y=self.dxy_h, boundary='trim').mean()

        ### Geometry
        
        print('*** Processing geometry')
        
        S, B, H, mask = Function(Q), Function(Q), Function(Q), Function(Q)
        kwargsLN = dict(x=xQ, y=yQ, method='linear')
        kwargsNN = dict(x=xQ, y=yQ, method='nearest')
        S.vector()[:]    = ds_b.surface.interp(**kwargsLN).to_numpy()
        B.vector()[:]    = ds_b.bed.interp(**kwargsLN).to_numpy()
        H.vector()[:]    = ds_b.thickness.interp(**kwargsLN).to_numpy()
        mask.vector()[:] = ds_b.mask.interp(**kwargsNN).to_numpy()

        ### Velocity
        
        print('*** Processing velocity')

        ds_ui = ds_u.interp(x=xQ2, y=yQ2, method='linear')
        ux_np = 1/self.velscale * getattr(ds_ui, self.uxname).to_numpy()
        uy_np = 1/self.velscale * getattr(ds_ui, self.uyname).to_numpy()
        I = np.argwhere(np.isnan(ux_np))
        if len(I) > 0:
            print('*** ERROR: NaNs found in velocity map where nodes are located...')
            if 1:
                print('...setting %i components to zero and silently proceeding.'%(len(I)))
                ux_np[I] = 1e-15
                uy_np[I] = 1e-15
            else:
                print('... exiting. I=', I.T, coordsQ2[:,I])
                sys.exit(1)
        ux, uy, u = Function(Q2), Function(Q2), Function(V)
        ux.vector()[:] = ux_np
        uy.vector()[:] = uy_np
        FunctionAssigner(V, [Q2,Q2]).assign(u, [ux,uy]) # set vector field components

        umag = Function(Q2)
        umag_np = np.sqrt(np.power(ux_np,2)+np.power(uy_np,2))
        umag.vector()[:] = umag_np
    
        D    = project(sym(grad(u)), TensorFunctionSpace(mesh, 'CG', 1))
        epsE = project(inner(D,D)+tr(D)**2, Q2)
        epsE.vector()[:] = np.sqrt(epsE.vector()[:]/2) # sqrt(D:D/2)

        if np.isnan(umag.vector()[:]).any(): raise ValueError('*** ERROR: umag contains NaNs')
        if (umag.vector()[:]<0).any():       raise ValueError('*** ERROR: umag contains <0 vals')

        if np.isnan(epsE.vector()[:]).any(): raise ValueError('*** ERROR: epsE contains NaNs')
        if (epsE.vector()[:]<0).any():       raise ValueError('*** ERROR: epsE contains <0 vals')

        ### Save FE fields

        print('*** Saving FE fields to inputs.h5')

        fout = '%s/inputs.h5'%(self.exname)
        f = HDF5File(MPI.comm_world, fout, "w")
        f.write(u, "/u")
        f.write(umag, "/umag")
        f.write(epsE, "/epsE")
        f.write(S, "/S")
        f.write(B, "/B")
        f.write(H, "/H")
        f.write(mask, "/mask")
        f.close()
        
        ### Save numpy fields

        print('*** Saving numpy fields to inputs.pkl')

        fout = '%s/inputs.pkl'%(self.exname)
        with open(fout, 'wb') as handle:
            data = [coords,cells] + [F.compute_vertex_values(mesh) for F in (ux,uy,umag,epsE, S,B,H,mask)]
            pickle.dump(data, handle)

    def get_mesh(self, scale=1):
        mesh = Mesh('%s/%s.xml'%(self.exname,self.meshname)) 
        boundaries = MeshFunction('size_t', mesh, '%s/%s_facet_region.xml'%(self.exname,self.meshname))
        Q  = FunctionSpace(mesh, 'CG', 1)
        Q2 = FunctionSpace(mesh, 'CG', 2) # this space is used when constructing the FEM velocity field (strain-rate tensor requires it to be DG2; Gibbs phenomenon?) 
        V  = VectorFunctionSpace(mesh, 'CG', 2)
        coords = scale*copy.deepcopy(mesh.coordinates().reshape((-1, 2)).T)
        cells = mesh.cells()
        return (mesh, boundaries, Q,Q2,V, coords,cells)
        
    def get_triang(self, coords, cells, scale=1):
        return tri.Triangulation(scale*coords[0], scale*coords[1], triangles=cells)
        
    def get_bmesh(self, scale=1):
        mesh, boundaries, *_ = self.get_mesh()
        boundarymesh = BoundaryMesh(mesh, 'exterior')
        bdim = boundarymesh.topology().dim()
        boundary_boundaries = MeshFunction('size_t', boundarymesh, bdim)
        boundary_boundaries.set_all(0)
        for i, facet in enumerate(entities(boundarymesh, bdim)):
            parent_meshentity = boundarymesh.entity_map(bdim)[i]
            parent_boundarynumber = boundaries.array()[parent_meshentity]
            boundary_boundaries.array()[i] = parent_boundarynumber
        bmesh_iso   = SubMesh(boundarymesh, boundary_boundaries, 1)
        bmesh_free  = SubMesh(boundarymesh, boundary_boundaries, 2)
        bcoords_iso  = copy.deepcopy(bmesh_iso.coordinates().reshape((-1, 2)).T)
        bcoords_free = copy.deepcopy(bmesh_free.coordinates().reshape((-1, 2)).T)
        return (bmesh_iso, bmesh_free, bcoords_iso, bcoords_free)
        
    def get_FE_inputs(self, Q=None, Q2=None, V=None):
        if Q is None: mesh,boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u,umag,epsE = Function(V), Function(Q2), Function(Q2)
        S,B,H,mask = Function(Q), Function(Q), Function(Q), Function(Q)
        f = HDF5File(MPI.comm_world,'%s/inputs.h5'%(self.exname),'r')
        f.read(u, '/u')
        f.read(umag, '/umag')
        f.read(epsE, '/epsE')
        f.read(S, '/S')
        f.read(B, '/B')
        f.read(H, '/H')
        f.read(mask, '/mask')
        f.close()
        return (u,umag,epsE, S,B,H,mask)
        
    def get_np_inputs(self):
        #print('*** Loading numpy inputs')
        with open('%s/inputs.pkl'%(self.exname), 'rb') as handle:
            return pickle.load(handle)
        
    """
    Fabric solver
    """

    def solve(self, FP, T=-15, altbc=False, verbose=True):
    
        ### Setup mesh and boundary conditions
    
        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u, umag, *_ = self.get_FE_inputs(Q=Q, Q2=Q2, V=V)
        
        nu_realspace = self.nu_realspace_DDRX if FP=='DDRX' else self.nu_realspace_LROT
        fab = IceFabric(mesh, boundaries, self.L, nu_realspace=nu_realspace, nu_multiplier=self.nu_multiplier, modelplane=self.modelplane, CAFFE_params=self.CAFFE_params)
        
        if altbc:
            x,y,z = np.eye(3) # Cartesian unit vectors
            sr = np.real(fab.sf.nlm_to_rnlm(fab.sf.nlm_ideal(y, 0, fab.L), fab.sf.get_rnlm_len()))
            fab.set_BCs([sr,], [0*sr,], [self.BOUNDARY_INFLOW,], domain=boundaries)
        else:
            fab.set_isotropic_BCs([self.BOUNDARY_INFLOW,], domain=boundaries) 

        ### Solve steady SSA fabric problem
        
        if FP not in ['LROT','DDRX']:
            raise ValueError('Invalid fabric process "%s"'%(FP))
        
        if verbose:
            if FP=='LROT': fabdynstr = 'LROT' 
            if FP=='DDRX': fabdynstr = 'LROT+DDRX' 
            print(tabulate([['Fabric dynamics', fabdynstr], 
                            ['Model plane', self.modelplane], 
                            ['L', self.L], 
                            ['nu_realspace', nu_realspace], 
                            ['nu_multiplier', self.nu_multiplier], 
                            ['T (deg. C)', T], 
                            ['Alternative BC?', altbc]
                            ], headers=['Parameter', 'Value'], tablefmt="rst"))
            print()
            
        print('*** Solving for steady SSA fabric field')

        if   FP == 'LROT': Gamma0 = None
        elif FP == 'DDRX': Gamma0 = [fab.Gamma0_Lilien23_lab(u, T+273.15) for T in np.linspace(-40,T,3)] # list of DDRX rate factors used to gradually approach solution 

        S = sym(grad(u)) # approximation strain rate tensor (D) as coaxial to stress tensor (S)
        fab.solvesteady(u, S, iota=+1, Gamma0=Gamma0, LROT_guess=True)
        
        ### Calculate fabric-derived quantities
        
        mi, lami = fab.eigenframe(*coords) # sorted such that last entry is the out-out-model-plane eigenpair
        E = fab.get_E_CAFFE(u).compute_vertex_values(mesh)    
        E[np.isnan(E)] = 1 # not sure why a few values sometimes are nan
            
        ### Save numpy solution

        suffix = '-altbc' if altbc else ''
        
        fout = '%s/solution-%s%s.pkl'%(self.exname, FP, suffix)
        print('*** Saving numpy solution to %s'%(fout))
        with open(fout, 'wb') as handle:
            data = (coords,cells, mi,lami,E)
            pickle.dump(data, handle)
         
        ### Save FE solution
        
        fout = '%s/solution-%s%s.h5'%(self.exname, FP, suffix)
        print('*** Saving FE solution to %s'%(fout))
        f = HDF5File(MPI.comm_world, fout, "w")
        f.write(fab.s,"/s")
        f.close()
        
    def get_np_solution(self, FP, suffix=''):
        fout = '%s/solution-%s%s.pkl'%(self.exname, FP, suffix)
        with open(fout, 'rb') as handle:
            #print('*** Loading solution %s'%(fout))
            return pickle.load(handle)
       
    def get_FE_solution(self, FP, suffix='', mesh=None, boundaries=None):
        if mesh is None: mesh, boundaries, *_ = self.get_mesh()
        fab = IceFabric(mesh, boundaries, self.L)
        s = Function(fab.S)
        fout = '%s/solution-%s%s.h5'%(self.exname, FP, suffix)
        f = HDF5File(MPI.comm_world,fout,"r")
        f.read(s,"/s")
        f.close()
        fab.s.assign(s)
        return fab
            
    """
    Plotting routines
    """
        
    def plot_inputs(self, lvls_u=None, lvls_D=None):

        (coords,cells, ux,uy,umag,epsE, S,B,H,mask) = self.get_np_inputs()
        triang = self.get_triang(coords, cells)

        fig, ax1, ax2 = self.setup_fig()
        
        ax = ax1
        F = np.log10(self.velscale*umag)
        lvls_u = np.linspace(np.percentile(F,5), np.percentile(F,95), 7) if lvls_u is None else lvls_u
        cs = ax.tricontourf(triang, F, levels=lvls_u, extend='both', locator=None, cmap='magma')
        hcb = plt.colorbar(cs, orientation='horizontal')
        hcb.set_label(r'$\log_{10} u$ (\SI{}{\metre\per\year})')
                
        ax = ax2
        F = self.velscale*epsE
        F = np.log10(self.velscale*epsE)
        lvls_D = np.linspace(np.percentile(F,1), np.percentile(F,90), 7) if lvls_D is None else lvls_D
        cs = ax.tricontourf(triang, F, levels=lvls_D, extend='both', locator=None, cmap='Spectral_r')
        hcb = plt.colorbar(cs, orientation='horizontal')
        hcb.set_label(r'$\log_{10} \dot{{\epsilon}}$ (\SI{}{\per\year})')

        for ax in (ax1,ax2):
            self.plot_mesh(ax, triang)
            self.plot_boundaries(ax)
            self.plot_floating(ax, triang, mask)

        fout = '%s/inputs.png'%(self.exname)
        print('*** Saving %s'%(fout))
        plt.savefig(fout, dpi=175, bbox_inches='tight')

    def plot_solution(self, FP, suffix='', Dlvls=np.linspace(0,1e-2,6)):

        (coords,cells, ux,uy,umag,epsE, S,B,H,mask) = self.get_np_inputs()
        (coords,cells, mi,lami,E_CAFFE) = self.get_np_solution(FP, suffix=suffix)
        triang = self.get_triang(coords, cells)
        
        fig, ax1, ax2 = self.setup_fig()
        
        ax = ax1
        dlam = abs(lami[:,0] - lami[:,1])
        cs = ax.tricontourf(triang, dlam, levels=np.arange(0,1.05,0.1), extend='neither', locator=None, cmap=plt.get_cmap('magma_r'))
        hcb = plt.colorbar(cs, orientation='horizontal')
        hcb.set_label(r'$\Delta \lambda$')
                
        ax = ax2
        lvlsE = np.arange(0, 3+1e-3, 0.25)
        divnorm = colors.TwoSlopeNorm(vmin=np.amin(lvlsE), vcenter=1, vmax=np.amax(lvlsE))
        kwargs_E = dict(levels=lvlsE, norm=divnorm, extend='max', cmap='PuOr_r')
        cs = ax.tricontourf(triang, E_CAFFE, **kwargs_E)
        hcb = plt.colorbar(cs, orientation='horizontal')
        hcb.set_label(r'$E$ (CAFFE)')
        
        for ax in (ax1,ax2):
            self.plot_mesh(ax, triang)
            self.plot_boundaries(ax)
            self.plot_floating(ax, triang, mask)

        fout = '%s/solution-%s%s.png'%(self.exname, FP, suffix)
        print('*** Saving %s'%(fout))
        plt.savefig(fout, dpi=175, bbox_inches='tight')
        
    def setup_fig(self):
        fig = plt.figure(figsize=(self.sx, self.sy))
        gs = gridspec.GridSpec(1,2)
        gs.update()
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        for ax in (ax1,ax2):
            self.plot_bg(ax)
            ax.axis('square')
            ax.set_xlabel(r'$x$ ($\SI{}{\metre}$)')
            ax.set_ylabel(r'$y$ ($\SI{}{\metre}$)')
            ax.set_xlim(self.xlims)
            ax.set_ylim(self.ylims)
        return fig, ax1, ax2
        
    def plot_mesh(self, ax, triang): 
        return ax.triplot(triang, **{'lw':0.075, 'color':'0.5', 'alpha':0.8})
    
    def plot_bg(self, ax):
        ax.add_patch(plt.Rectangle((self.x0,self.y0), np.diff(self.xlims)[0], np.diff(self.ylims)[0], color='0.85'))
    
    def plot_boundaries(self, ax, scale=1, lw=2, zorder=20, clip_on=False, 
                              c_inflow='c', c_free='#f0027f', ls_inflow='-', ls_free='--', lbl_inflow='Isotropic', lbl_free='Free', **kwargs):

        (bmesh_iso, bmesh_free, bcoords_iso, bcoords_free) = self.get_bmesh()
        xy_iso  = bcoords_iso.T
        xy_free = bcoords_free.T
        tour_iso,  total_distance = self.tsp(xy_iso)
        tour_free, total_distance = self.tsp(xy_free)
        x_iso  = np.array([xy_iso[i][0]  for i in tour_iso])
        y_iso  = np.array([xy_iso[i][1]  for i in tour_iso])
        x_free = np.array([xy_free[i][0] for i in tour_free])
        y_free = np.array([xy_free[i][1] for i in tour_free])
        h_iso  = ax.plot(x_iso*scale,  y_iso*scale,  ls=ls_inflow,  c=c_inflow,  lw=lw, label=lbl_inflow, zorder=zorder, clip_on=clip_on, **kwargs)
        h_free = ax.plot(x_free*scale, y_free*scale, ls=ls_free,    c=c_free,    lw=lw, label=lbl_free,   zorder=zorder, clip_on=clip_on, **kwargs)
        return (h_iso, h_free)
    
    def plot_floating(self, ax, triang, mask, scale=1, c='limegreen', lw=1.8):
        ax.tricontour(triang, mask==3, [0.5, 1.5], colors=[c,], linewidths=lw)
        h_float = ax.plot([0,0], [1e10,1e10], c=c, lw=lw, label='Floating') # fake entry for legend
        return h_float

    """
    AUX
    """
        
#    def plot_stream(ax, Xv, Yv, ux,uy,meshpts):
#    
#        velcutoff = 0 # don't draw streamlines when velocity is lower than this threshold (m/yr)
#        umag = np.sqrt(np.power(ux,2)+np.power(uy,2))
#        ux[umag < velcutoff] = np.nan # no stream lines over low-velocity areas to avoid noisy plot
#        X2, Y2 = np.meshgrid(Xv, Yv, indexing='xy')
#        UX2 = griddata(meshpts, ux.flatten(), (X2, Y2), method='linear', fill_value=np.nan)
#        UY2 = griddata(meshpts, uy.flatten(), (X2, Y2), method='linear', fill_value=np.nan)

#        x0, y0 = 195, -1500
#        x1, y1 = 340, -1550
#        k = (y1-y0)/(x1-x0)
#        sp = [(x,y0+k*(x-x0)) for x in np.linspace(x0, x1, 12)]
#        kwargs_streamplot = dict(density=100, linewidth=1, color='0.9', broken_streamlines=False, start_points=sp)

#        ax.streamplot(X2, Y2, UX2, UY2, **kwargs_streamplot)

#        conf.plot_velocities(ax, umag, triang)
#        conf.plot_mesh(ax, triang)
#        conf.plot_stream(ax, sc*Xv, sc*Yv, ux,uy,meshpts)
#        if conf.debug: conf.plot_boundaries(ax, MASK, X*sc, Y*sc, legend=True)


    """
    TSP for determining ordered list of boundary coordinates 
    Code from https://www.askpython.com/python/examples/travelling-salesman-problem-python
    """
    
    def distance(self, city1, city2):
      # Replace this with your distance calculation function (e.g., Euclidean distance)
      x1, y1 = city1
      x2, y2 = city2
      return ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
     
    def tsp(self, cities):
      visited = [False] * len(cities)
      current_city = 0
     
      tour = []
      total_distance = 0
     
      for _ in range(len(cities)):
        visited[current_city] = True
        tour.append(current_city)
     
        next_city = None
        min_distance = float('inf')  # Represents positive infinity
     
        for i in range(len(cities)):
          if visited[i]:
            continue
     
          d = self.distance(cities[current_city], cities[i])
          if d < min_distance:
            min_distance = d
            next_city = i
     
        current_city = next_city
        total_distance += min_distance
     
      return tour, total_distance
      
