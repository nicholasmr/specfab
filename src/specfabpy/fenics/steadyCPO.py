#!/usr/bin/python3
# Nicholas Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Steady SSA CPO solver for ice
"""

import os, sys, copy, code # code.interact(local=locals())
import numpy as np
import xarray, pickle, pyproj # pyproj needed!
from tabulate import tabulate

#from dolfin import *
from .ice import * # IceFabric

import matplotlib.tri as tri
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

class steadyCPO():

    uxname, uyname = 'VX', 'VY'
    xuname, yuname = 'x', 'y'

    CAFFE_params = (0.1, 10) # Emin, Emax of CAFFE

    ms2myr = 3.17098e+8 # m/s to m/yr
    m2km = 1e-3 # m to km

    modelplane = 'xz' # xz is much faster than xy (renders state vector real-valued)

    bc_isotropic  = 0
    bc_zsinglemax = 1

    c_floating = 'limegreen'
    
    kw_leg = dict(ncol=3, loc=1, bbox_to_anchor=(1,1.12), handlelength=1.2, fancybox=False, frameon=False)

    def __init__(self, domain):
    
        self.exname = domain['name'] # experiment name
        self.fvel   = domain['fmeasures']
        self.fbed   = domain['fbedmachine']
        self.fmesh  = domain['fgeo'][:-4]
        self.dxy_u  = domain['subsample_u']
        self.dxy_h  = domain['subsample_h']
        
        ### Determine bounding box of domain for interpolation

        content = open(domain['fgeo']).readlines()
        p0 = (self._getnum(content[2]), self._getnum(content[3]))
        p1 = (self._getnum(content[4]), self._getnum(content[5]))
        #print(p0,p1)
        # ...unpack
        self.x0, self.y0 = p0
        self.x1, self.y1 = p1
        self.dx = 0.05 * abs(self.x1-self.x0)
        self.dy = 0.05 * abs(self.y1-self.y0)
        
        self.mapscale = self.m2km # x and y axis scale
        
    def _getnum(self, s):
        part1, part2 = s.split("=", 1) # Split the string into two parts based on the start delimiter
        result = part2.split(";", 1)[0] # Split the second part based on the end delimiter
        return float(result)
        
    def preprocess(self):

        """
        Slice and subsample maps of velocity and geometry.
        Results are saved to .h5 and .pkl files.
        """
    
        print('*** Preprocessing started')
        print(tabulate([['name', self.exname], ['fmeasures', self.fvel], ['fbedmachine', self.fbed], \
                        ['subsample_u', self.dxy_u], ['subsample_h', self.dxy_h]], \
                        headers=['Parameter', 'Value'], tablefmt="rst"))

        ### Make mesh

        print('*** Generating finite element mesh')
        os.system('gmsh -2 -format msh2 %s.geo'%(self.fmesh))
        os.system('dolfin-convert %s.msh %s.xml'%(self.fmesh, self.fmesh))

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

        print('*** Coarsening input fields')

        ds_vel0 = xarray.open_mfdataset(self.fvel)
        ds_bed0 = xarray.open_mfdataset(self.fbed)
        xrng = slice(self.x0-self.dx, self.x1+self.dx)
        yrng = slice(self.y1+self.dy, self.y0-self.dy)
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
        ux_np = 1/self.ms2myr * getattr(ds_ui, self.uxname).to_numpy()
        uy_np = 1/self.ms2myr * getattr(ds_ui, self.uyname).to_numpy()
        I = np.argwhere(np.isnan(ux_np))
        if len(I) > 0:
            print('*** Warning: NaNs found in velocity map where nodes are located... please make sure the model domain does not extend outside that of the velocity product')
            if 1:
                print('...setting %i components to zero and proceeding silently.'%(len(I)))
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

        if np.isnan(umag.vector()[:]).any(): raise ValueError('*** Warning: umag contains NaNs')
        if (umag.vector()[:]<0).any():       raise ValueError('*** Warning: umag contains values <0 ')

        if np.isnan(epsE.vector()[:]).any(): raise ValueError('*** Warning: epsE contains NaNs')
        if (epsE.vector()[:]<0).any():       raise ValueError('*** Warning: epsE contains values <0 ')

        ### Save FEM fields

        print('*** Saving FEM fields to %s-inp.h5'%(self.exname))

        fout = '%s-inp.h5'%(self.exname)
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

        print('*** Saving numpy fields to %s-inp.pkl'%(self.exname))

        fout = '%s-inp.pkl'%(self.exname)
        with open(fout, 'wb') as handle:
            data = [coords,cells] + [F.compute_vertex_values(mesh) for F in (ux,uy,umag,epsE, S,B,H,mask)]
            pickle.dump(data, handle)

    def get_mesh(self, mapscale=1):
    
        mesh = Mesh('%s.xml'%(self.fmesh)) 
        boundaries = MeshFunction('size_t', mesh, '%s_facet_region.xml'%(self.fmesh))
        Q  = FunctionSpace(mesh, 'CG', 1)
        Q2 = FunctionSpace(mesh, 'CG', 2) # this space is used when constructing the FEM velocity field (strain-rate tensor requires it to be DG2; Gibbs phenomenon?) 
        V  = VectorFunctionSpace(mesh, 'CG', 2)
        coords = mapscale*copy.deepcopy(mesh.coordinates().reshape((-1, 2)).T)
        cells = mesh.cells()
        return (mesh, boundaries, Q,Q2,V, coords,cells)
        
    def triang(self, coords, cells, mapscale=1):
    
        return tri.Triangulation(mapscale*coords[0], mapscale*coords[1], triangles=cells)
        
    def feminputs(self, Q=None, Q2=None, V=None):
    
        if Q is None: mesh,boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u,umag,epsE = Function(V), Function(Q2), Function(Q2)
        S,B,H,mask = Function(Q), Function(Q), Function(Q), Function(Q)
        f = HDF5File(MPI.comm_world,'%s-inp.h5'%(self.exname),'r')
        f.read(u, '/u')
        f.read(umag, '/umag')
        f.read(epsE, '/epsE')
        f.read(S, '/S')
        f.read(B, '/B')
        f.read(H, '/H')
        f.read(mask, '/mask')
        f.close()
        return (u,umag,epsE, S,B,H,mask)
        
    def npinputs(self):
    
        #print('*** Loading numpy inputs')
        with open('%s-inp.pkl'%(self.exname), 'rb') as handle:
            return pickle.load(handle)
        
    """
    Solver
    """

    def solve(self, problem, numerics):
    
        isDDRX = len(np.shape(problem['T']))>0
        print(tabulate([['PROBLEM', 'LROT+ADVEC' if not isDDRX else 'LROT+DDRX+ADVEC'], 
                        ['T',          ','.join(['%.1f'%(_) for _ in problem['T']]) if isDDRX else 'None'], 
                        ['A,Q (DDRX)', ','.join(['%.2e'%(_) for _ in problem['AQ_DDRX']]) if 'AQ_DDRX' in problem.keys() else 'None'], 
                        ['L',          numerics['L']], 
                        ['nu_orimul',  '%.2e'%(numerics['nu_orimul'])], 
                        ['nu_real',    '%.2e'%(numerics['nu_real'])], 
                        ['nu_realmul', '1' if not isDDRX else ','.join(['%.1f'%(_) for _ in numerics['nu_realmul']])], 
                        ['modelplane', self.modelplane], 
                        ], headers=['Parameter', 'Value'], tablefmt="rst"))
    
        ### Mesh

        print('*** Loading mesh and input fields')
    
        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u, umag, *_ = self.feminputs(Q=Q, Q2=Q2, V=V)        
        
        ### Initialize solver class
        
        fab = IceFabric(mesh, boundaries, numerics['L'], nu_realspace=numerics['nu_real'], nu_multiplier=numerics['nu_orimul'], \
                            modelplane=self.modelplane, CAFFE_params=self.CAFFE_params)
        
        ### Boundary conditions  
        
        nlm_list = [ [1/np.sqrt(4*np.pi)] + [0]*(fab.nlm_len_full-1), fab.sf.nlm_ideal([0,1,0], 0, numerics['L'])] # isotropic or z-SMAX (note y is vertical since modelplane is xz)
        bc_vals  = [np.real(fab.sf.nlm_to_rnlm(nlm_list[bc[1]], fab.nlm_len)) for bc in problem['bcs']]
        bc_ids   = [bc[0] for bc in problem['bcs']]
        fab.set_BCs(bc_vals, [0*val for val in bc_vals], bc_ids, domain=boundaries) # (real, imag) parts on bc_ids
            
        ### Solve steady SSA problem
        
        if not isDDRX: Gamma0 = None
        else:          Gamma0 = [fab.Gamma0(u, _+273.15, *problem['AQ_DDRX']) for _ in problem['T']] # list of DDRX rate factors used to gradually approach solution 

        S = sym(grad(u)) # approximation strain rate tensor (D) as coaxial to stress tensor (S)
        fab.solvesteady(u, S, iota=+1, Gamma0=Gamma0, nu_realspacemul=numerics['nu_realmul'], LROT_guess=True)
        
        ### Calculate fabric-derived quantities
        
        mi, lami = fab.eigenframe(*coords) # sorted such that last entry is the out-out-model-plane eigenpair
        E = fab.get_E_CAFFE(u).compute_vertex_values(mesh)    
        E[np.isnan(E)] = 1 # not sure why a few values sometimes are nan
            
        ### Save numpy solution

        fout = '%s-sol-%s.pkl'%(self.exname, problem['name'])
        print('*** Saving numpy solution to %s'%(fout))
        with open(fout, 'wb') as handle:
            data = (coords,cells, mi,lami,E)
            pickle.dump(data, handle)
         
        ### Save FEM solution
        
        fout = '%s-sol-%s.h5'%(self.exname, problem['name'])
        print('*** Saving FEM solution to %s'%(fout))
        f = HDF5File(MPI.comm_world, fout, "w")
        f.write(fab.s,"/s")
        f.close()
        
    def npsolution(self, probname):
    
        fout = '%s-sol-%s.pkl'%(self.exname, probname)
        with open(fout, 'rb') as handle:
            return pickle.load(handle)
       
    def femsolution(self, probname, L, mesh=None, boundaries=None):
        
        if mesh is None: mesh, boundaries, *_ = self.get_mesh()
        self.fab = IceFabric(mesh, boundaries, L)
        s = Function(self.fab.S)
        fout = '%s-sol-%s.h5'%(self.exname, probname)
        f = HDF5File(MPI.comm_world,fout,"r")
        f.read(s,"/s")
        f.close()
        self.fab.s.assign(s)
        return self.fab
            
    """
    AUX
    """

    def set_solution(self, probname):
    
        self.coords,self.cells, self.mi,self.lami,self.E_CAFFE = self.npsolution(probname)

    def set_inputs(self):
    
        self.coords,self.cells, self.ux,self.uy,self.umag,self.epsE, self.S,self.B,self.H,self.mask = self.npinputs()
        self.triang = self.triang(self.coords, self.cells, mapscale=self.mapscale)

    def bmesh(self, mapscale=1):
    
        mesh, boundaries, *_ = self.get_mesh()
        boundarymesh = BoundaryMesh(mesh, 'exterior')
        bdim = boundarymesh.topology().dim()
        boundary_boundaries = MeshFunction('size_t', boundarymesh, bdim)
        boundary_boundaries.set_all(0)
        for i, facet in enumerate(entities(boundarymesh, bdim)):
            parent_meshentity = boundarymesh.entity_map(bdim)[i]
            parent_boundarynumber = boundaries.array()[parent_meshentity]
            boundary_boundaries.array()[i] = parent_boundarynumber
        bmeshes = [SubMesh(boundarymesh, boundary_boundaries, bcid) for bcid in [1,2]] # assume only two possible kinds; isotropic/smax and free
        coords = [copy.deepcopy(bmesh.coordinates().reshape((-1, 2)).T) for bmesh in bmeshes]
        return (coords, bmeshes)
        
    """
    Plotting
    """

    def plot_inputs(self, figsize=(5,5), boundaries=False, mesh=True, kw_vel={}, kw_epsE={}):
    
        self.set_inputs()
        
        fig, ax = self.newfig(figsize=figsize)
        legh, legt = self.setupaxis(ax, boundaries=boundaries, mesh=mesh)
        if boundaries: ax.legend(legh, legt, **self.kw_leg)
        self.plot_velocities(ax, **kw_vel)
        self.savefig(fig, '%s-inp-vel.png'%(self.exname))
        
        fig, ax = self.newfig(figsize=figsize)
        legh, legt = self.setupaxis(ax, boundaries=boundaries, mesh=mesh)
        if boundaries: ax.legend(legh, legt, **self.kw_leg)
        self.plot_strainratemag(ax, **kw_epsE)
        self.savefig(fig, '%s-inp-epsE.png'%(self.exname))
      
    def plot_results(self, problem, figsize=(5,5), boundaries=True, mesh=False, kw_E={}, kw_dlam={}, kw_lamz={}):
    
        self.set_inputs()
        self.set_solution(problem['name'])
        kw_setupaxis = dict(boundaries=boundaries, mesh=mesh)
        
        fig, ax = self.newfig(figsize=figsize)
        legh, legt = self.setupaxis(ax, **kw_setupaxis)
        ax.legend(legh, legt, **self.kw_leg)
        self.plot_E_CAFFE(ax, **kw_E)
        self.savefig(fig, '%s-sol-%s-E.png'%(self.exname, problem['name']))
            
        fig, ax = self.newfig(figsize=figsize)
        legh, legt = self.setupaxis(ax, **kw_setupaxis)
        ax.legend(legh, legt, **self.kw_leg)
        self.plot_dlam(ax, **kw_dlam)
        self.savefig(fig, '%s-sol-%s-dlam.png'%(self.exname, problem['name']))
        
        fig, ax = self.newfig(figsize=figsize)
        legh, legt = self.setupaxis(ax, **kw_setupaxis)
        ax.legend(legh, legt, **self.kw_leg)
        self.plot_lamz(ax, **kw_lamz)
        self.savefig(fig, '%s-sol-%s-lamz.png'%(self.exname, problem['name']))

    def plot_generic(self, ax, F, kw_tcf=dict(), kw_cb=dict(), kw_cax=dict()):

        if ('norm' in kw_tcf) and (kw_tcf['norm'] == 'log'):
            kw_tcf['norm'] = colors.LogNorm(vmin=kw_tcf['levels'][0], vmax=kw_tcf['levels'][-1])    
        if ('norm' in kw_tcf) and (kw_tcf['norm'] == 'center'): 
            kw_tcf['norm'] = colors.TwoSlopeNorm(vmin=kw_tcf['levels'][0], vcenter=kw_tcf['vcenter'], vmax=kw_tcf['levels'][-1])
            kw_tcf.pop('vcenter')
        cs = ax.tricontourf(self.triang, F, **kw_tcf)
        hcb = plt.colorbar(cs, cax=self.newcax(ax, **kw_cax), **kw_cb)
        return (cs, hcb)
        
    def plot_velocities(self, ax, kw_cb=dict(label=r'$u$ (m/yr)'), kw_cax=dict(), \
                            kw_tcf=dict(cmap='inferno', levels=np.logspace(0.5, 3.5, 13), norm='log', extend='both')):
                            
        return self.plot_generic(ax, self.ms2myr*self.umag, kw_tcf=kw_tcf, kw_cb=kw_cb, kw_cax=kw_cax)

    def plot_strainratemag(self, ax, kw_cb=dict(label=r'$\dot{\epsilon}_{e}$ (1/yr)'), kw_cax=dict(), \
                            kw_tcf=dict(cmap='viridis', levels=np.arange(0, 50+.01, 5), extend='max')):
                            
        return self.plot_generic(ax, 1e3*self.ms2myr*self.epsE, kw_tcf=kw_tcf, kw_cb=kw_cb, kw_cax=kw_cax)

    def plot_dlam(self, ax,  kw_cb=dict(label=r'$\Delta\lambda$'), kw_cax=dict(), \
                        kw_tcf=dict(cmap='Blues', levels=np.arange(0, 0.8+.01, 0.1), extend='max'), \
                        quiver=False, quiverkey=(0.1, 0.05, 3, r'${\bf m}_1$') ): 

        dlam = abs(self.lami[:,0] - self.lami[:,1]) # eigenvalues 1 and 2 are the largest and smallest in-model-plane eigenvalues
        returnme = self.plot_generic(ax, dlam, kw_tcf=kw_tcf, kw_cb=kw_cb, kw_cax=kw_cax)

        # Quiver principal horizontal eigenvector?
        if quiver:
            meshpts = (coords[0,:], coords[1,:])
            xv, yv = np.linspace(scpo.x0, scpo.x1, 15)[1:-1], np.linspace(scpo.y0, scpo.y1, 15)[1:-1]
            x, y = np.meshgrid(xv, yv, indexing='xy')
            m1 = mi[:,0,:] # principal horizontal eigenvector
            m1x = griddata(meshpts, m1[:,0].flatten(), (x, y), method='linear', fill_value=np.nan)
            m1y = griddata(meshpts, m1[:,2].flatten(), (x, y), method='linear', fill_value=np.nan) # y coordinate is index 2 (z coordinate) since problem is in xz plane
            renorm = np.sqrt(m1x**2+m1y**2)
            m1x, m1y = np.divide(m1x, renorm), np.divide(m1y, renorm)
            hq = ax.quiver(mapscale*x, mapscale*y, +m1x, +m1y, color='tab:red', scale=40)
            hq = ax.quiver(mapscale*x, mapscale*y, -m1x, -m1y, color='tab:red', scale=40)
            ax.quiverkey(hq, *quiverkey, labelpos='E')

        return returnme

    def plot_lamz(self, ax, kw_cb=dict(label=r'$\lambda_z$'), kw_cax=dict(), \
                        kw_tcf=dict(cmap='RdPu', levels=np.arange(0, 0.8+.01, 0.1), extend='max')): 
                        
        return self.plot_generic(ax, self.lami[:,2], kw_tcf=kw_tcf, kw_cb=kw_cb, kw_cax=kw_cax) # eigenvalue 3 is the out-of-model-plane (z) eigenvalue

    def plot_E_CAFFE(self, ax, kw_cb=dict(label=r'$E$'), kw_cax=dict(), \
                            kw_tcf=dict(cmap='PuOr_r', levels=np.logspace(-1, 1, 17), norm='log', extend='both')): 
                            
        return self.plot_generic(ax, self.E_CAFFE, kw_tcf=kw_tcf, kw_cb=kw_cb, kw_cax=kw_cax)
        
    def plot_shearfrac(self, ax, kw_cb=dict(label=r'{\fontsize{10}{10}\selectfont $\leftarrow$ stretching} \;\; $\gamma$\;\; {\fontsize{10}{10}\selectfont shearing $\rightarrow$}'), kw_cax=dict(), \
                            kw_tcf=dict(cmap='Spectral_r', levels=np.arange(0, 1+.01, 0.1), extend='neither')): 

        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u, *_ = self.feminputs(Q=Q, Q2=Q2, V=V)
        u.set_allow_extrapolation(True)
        shearfrac_df, *_ = self.fab.enhancementfactor.shearfrac_SSA(u)
        self.shearfrac   = shearfrac_df.compute_vertex_values(mesh)
        return self.plot_generic(ax, self.shearfrac, kw_tcf=kw_tcf, kw_cb=kw_cb, kw_cax=kw_cax)
    
    def plot_chi(self, ax, kw_cb=dict(label=r'{\fontsize{10}{10}\selectfont $\leftarrow$ incompatible} \;\; $\chi$\;\; {\fontsize{10}{10}\selectfont compatible $\rightarrow$}'), kw_cax=dict(), \
                            kw_tcf=dict(cmap='PiYG', levels=np.arange(0,1+.01,0.1), extend='neither')): 

        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u, *_ = self.feminputs(Q=Q, Q2=Q2, V=V)
        u.set_allow_extrapolation(True)
        chi_df   = self.fab.enhancementfactor.chi(u, self.fab.s)
        self.chi = chi_df.compute_vertex_values(mesh)
        return self.plot_generic(ax, self.chi, kw_tcf=kw_tcf, kw_cb=kw_cb, kw_cax=kw_cax)
        
    def newfig(self, figsize=(5,5)):
    
        fig = plt.figure(figsize=figsize)
        ax = plt.subplot(111)
        return (fig, ax)
        
    def savefig(self, fig, fname, dpi=150, pad_inches=0.1, bbox_inches='tight'):
    
        fig.savefig(fname, dpi=dpi, pad_inches=pad_inches, bbox_inches=bbox_inches)
        
    def setupaxis(self, ax, boundaries=False, floating=True, mesh=False, bgcolor='0.85', showyaxis=True, \
             xlims=None, ylims=None, xticks_major=None, xticks_minor=None, yticks_major=None, yticks_minor=None):

        legh, legt = [], []
        
        if bgcolor is not None: 
            x0, x1 = self.mapscale*(self.x0-self.dx), self.mapscale*(self.x1+self.dx)
            y0, y1 = self.mapscale*(self.y0-self.dy), self.mapscale*(self.y1+self.dy)
            ax.add_patch(plt.Rectangle((x0,y0), x1-x0, y1-y0, color=bgcolor))
            
        if mesh: 
            ax.triplot(self.triang, lw=0.075, color='0.5', alpha=0.8, zorder=10)
            
        if boundaries:
            coords, *_ = self.bmesh()
            colors = ['0.3', 'deeppink'] # , 'yellow', 'aquamarine']
            markers = ['s',]*len(colors)
            for ii, (xb, yb) in enumerate(coords):
                ax.scatter(xb*self.mapscale, yb*self.mapscale, c=colors[ii], marker=markers[ii], s=3, zorder=12, clip_on=False)
                legh.append(Line2D([0], [0], color=colors[ii], lw=2))
                legt.append('Isotropic' if ii == 0 else 'Free')
                
        if floating: 
            ax.tricontour(self.triang, self.mask==3, [0.5, 1.5], colors=[self.c_floating,], linewidths=2, zorder=11)
            legh.append(Line2D([0], [0], color=self.c_floating, lw=2))
            legt.append('Floating')
            
        ax.axis('square')
        ax.set_xlabel(r'$x$ (km)')

        if xticks_major is not None: ax.set_xticks(xticks_major)
        if xticks_minor is not None: ax.set_xticks(xticks_minor, minor=True)
            
        if showyaxis:
            if yticks_major is not None: ax.set_yticks(yticks_major)
            if yticks_minor is not None: ax.set_yticks(yticks_minor, minor=True)
            ax.set_ylabel(r'$y$ (km)')
        else:
            ax.tick_params('y', labelleft=False)

        ax.set_xlim([self.mapscale*self.x0, self.mapscale*self.x1] if xlims is None else xlims)
        ax.set_ylim([self.mapscale*self.y0, self.mapscale*self.y1] if ylims is None else ylims)
        
        return (legh, legt)
    
    def newcax(self, ax, loc="right", size="4%", pad=0.13, **kwargs): 
        return make_axes_locatable(ax).append_axes(loc, size=size, pad=pad, **kwargs)

