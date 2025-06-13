#!/usr/bin/python3
# Nicholas Rathmann and Daniel Shapero, 2024

r"""
Test firedrake interface for specfab.

Assumes a time-constant, non-uniform shear flow in vertical cross-section (xz) domain.
"""

import copy, sys, time, code # code.interact(local=locals())
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib.tri as tri
from matplotlib import rc
rc('font',**{'family':'serif', 'sans-serif':['Times'], 'size':15})

import firedrake as fd
from specfabpy.firedrake.ice import IceFabric
from specfabpy import plotting as sfplt

"""
Fabric problem setup
"""

### Numerics and regularization

Nt = 20 # number of time steps to take
L  = 8  # spectral truncation *** set L=8 or L=10 unless debugging ***

kw_num = dict(nu_multiplier=1, nu_realspace=1e-3, modelplane='xz') # nu_realspace must be adjusted when changing mesh/resolution (trial-and-error so that solution is stable and smooth)

DEBUG__STEADYSTATE = False # Do not solve time-dependent problem but show steady state solution

### Fabric dynamics

iota = +1  # deck-of-cards behaviour for lattice rotation
Tice = -15 # Ice temperature (deg. C) for DDRX rate factor

ENABLE_DDRX = True

### Viscous anisotropy homogenization parameters

alpha     = 0.455     # Taylor--Sachs homogenization weight
Eij_grain = (1, 1e3)  # (Ecc, Eca) grain enhancements
n_grain   = 1         # grain power-law exponent (only n_grain=1 supported)
E_CAFFE   = (0.1, 10, 1) # (Emin, Emax, n_grain) of CAFFE

kw_VA = dict(alpha=alpha, Eij_grain=Eij_grain, n_grain=n_grain, E_CAFFE=E_CAFFE)

"""
Setup firedrake fabric class
"""

### Mesh and function spaces

nx = ny = 16
mesh = fd.UnitSquareMesh(nx, ny, diagonal='right') #, diagonal='crossed')
boundaries = (1,2,3,4) 
x = fd.SpatialCoordinate(mesh)
V = fd.VectorFunctionSpace(mesh, "CG", 1)
T = fd.TensorFunctionSpace(mesh, "CG", 2)
Q = fd.FunctionSpace(mesh, "CG", 1) # for projecting scalar fabric measures

### Velocity field

u0, H = 1, 1
expr = fd.as_vector(( u0*(x[1]/H)**2, 0 )) # non-uniform horizontal shear
u = fd.Function(V).interpolate(expr)
tau = fd.project(fd.sym(fd.grad(u)), T) # assume driving stress is coaxial to strain-rate (in two-way coupling this should be modelled tau) 

h_min = 1/nx
v_max = abs(u.vector()[:]).max()
dt_CFL = 0.5*h_min/v_max
dt = 2*dt_CFL # more aggresive time-stepping than CFL

### Initialize fabric module

fabric = IceFabric(mesh, boundaries, L, **kw_num, **kw_VA) # initializes as isotropic fabric field
fabric.set_isotropic_BCs((1,)) # isotropic ice incoming from left-hand boundary, remaining boundaries are free (no fabric fluxes)

Gamma0 = None
if ENABLE_DDRX:
    Gamma0      = fabric.Gamma0_Lilien23_EDC(u,Tice+273.15)
    Gamma0_list = [fabric.Gamma0_Lilien23_EDC(u,_+273.15) for _ in np.linspace(-40,Tice,3)] # list of DDRX rate factors used to gradually approach solution 

"""
Solve for steady state 
"""

fabric.solvesteady(u, tau, iota=iota, Gamma0=Gamma0_list, LROT_guess=True)
pfJ_steady = fabric.get_pfJ().copy(deepcopy=True)

"""
Time evolution
"""

if not DEBUG__STEADYSTATE: 
    fabric.initialize() # reset to isotropic after solving for steady state

nn = 0
t  = 0.0

while nn < Nt:

    nn += 1
    t  += dt
    
    if not DEBUG__STEADYSTATE:
        print("*** Step %i :: dt=%.2e, t=%.2e" % (nn, dt, t))
        fabric.evolve(u, tau, dt, iota=iota, Gamma0=Gamma0) # automatically updates derived properties (Eij, pfJ, etc.)
    
    ### Plot results
        
    if nn==1 or nn%1==0:
    
        fname = 'simpleshear-%s-%i.png'%('LROT' if Gamma0 is None else 'DDRX', nn) # %02d
        print('... plotting model state: %s'%(fname))

        ### Setup figure
        
        figscale = 1.2
        figsize = (20*figscale, 9*figscale)
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 6, wspace=0.45, hspace=0.1, left=0.04, right=0.97, top=0.98, bottom=0.27)
        axr1 = [fig.add_subplot(gs[0,ii]) for ii in range(6)]
        axr2 = [fig.add_subplot(gs[1,ii]) for ii in range(6)]
        axes = np.concatenate((axr1,axr2))
        
        kwargs_cb = dict(aspect=22, fraction=0.2, orientation='horizontal')

        ### Plot J index

        ax = axr1[0]
        h = fd.pyplot.tricontourf(fd.project(fabric.pfJ,Q), axes=ax, levels=np.arange(1, 3+1e-3, 0.2), extend='max', cmap='YlGnBu')
        cbar = plt.colorbar(h, ax=ax, **kwargs_cb)
        cbar.ax.set_xlabel(r'$J$ index')
        h = fd.pyplot.quiver(u, axes=ax, cmap='Reds', width=0.0075)
        
        ax = axr1[1]
        h = fd.pyplot.tricontourf(fd.project(pfJ_steady,Q), axes=ax, levels=np.arange(1, 3+1e-3, 0.2), extend='max', cmap='YlGnBu')
        cbar = plt.colorbar(h, ax=ax, **kwargs_cb)
        cbar.ax.set_xlabel(r'$J$ index (steady-state)')
        h = fd.pyplot.quiver(u, axes=ax, cmap='Reds', width=0.0075)

        ax = axr1[2]
        lvls_E = np.arange(0.6, 2+1e-3, 0.1)
        divnorm_E = colors.TwoSlopeNorm(vmin=np.amin(lvls_E), vcenter=1, vmax=np.amax(lvls_E))
        kwargs_E = dict(levels=lvls_E, norm=divnorm_E, extend='both', cmap='PuOr_r')
        h = fd.pyplot.tricontourf(fd.project(fabric.E_CAFFE,Q), axes=ax, **kwargs_E)
        cbar = plt.colorbar(h, ax=ax, **kwargs_cb)
        cbar.ax.set_xlabel(r'$E$ (CAFFE)')
        
        c = ['Reds', 'Greens', 'Blues']
        for ii, ax in enumerate(axr1[3:]):
            h = fd.pyplot.tricontourf(fd.project(fabric.lami[ii],Q), axes=ax, levels=np.arange(0.0, 1+1e-3, 0.1), extend='neither', cmap=c[ii])
            cbar = plt.colorbar(h, ax=ax, **kwargs_cb)
            cbar.ax.set_xlabel(r'$\lambda_%i$'%(ii+1))

        idx = ['11','22','33','23','13','12'] # Voigt ordering
        for ii, ax in enumerate(axr2[:]):
            h = fd.pyplot.tricontourf(fd.project(fabric.Eij[ii],Q), axes=ax, **kwargs_E)
            cbar = plt.colorbar(h, ax=ax, **kwargs_cb)
            cbar.ax.set_xlabel(r'$E_{%s}$'%(idx[ii]))
       
        # plot grids 
        for ax in axes:
            ax.set_aspect('equal')
            ax.set_xlabel('$x$')
            ax.set_ylabel('$z$')
            ax.set_xlim([0,1])
            ax.set_ylim([0,1])
            fd.pyplot.triplot(mesh, axes=ax, interior_kw=dict(lw=0.1), boundary_kw=dict(lw=5, clip_on=False))

        ### ODF insets
        
        def plot_ODF(p, pax, mrk, W=0.15):
            
            geo, prj = sfplt.getprojection(rotation=-110, inclination=50)
            axin = plt.axes([*pax, W,W], projection=prj) 
            axin.set_global()
            sfplt.plotODF(fabric.get_nlm(*p), fabric.lm, axin, cmap='Greys', cblabel='ODF', lvlset=(np.linspace(0.0, 0.3, 6), lambda x,p:'%.1f'%x), showcb=True)
            sfplt.plotcoordaxes(axin, geo, color='k')
            mi, lami = fabric.eigenframe(*p)
            sfplt.plotmi(axin, mi, geo, ms=15, colors=(sfplt.c_dred,sfplt.c_dgreen,sfplt.c_dblue))
            axin.set_title(r'@ "%s"'%(mrk))
            for ax in axes: points, = ax.plot(*p, mrk, markersize=12, markeredgewidth=1.1, markeredgecolor='k', markerfacecolor='w')

        # Plot for selected locations
        plot_ODF((0.1,0.1), (0.3,0.07), '^')
        plot_ODF((0.5,0.4), (0.4,0.07), 'X')
        plot_ODF((0.9,0.6), (0.5,0.07), 's')
            
        axr2[4].text(-0.6, -1.2, '${\\bf m}_1$ is red\n${\\bf m}_2$ is green\n${\\bf m}_3$ is blue', ha='left', ma='left')
            
        ### Save plot
        
        plt.savefig(fname, dpi=100)
        plt.close()
        print('... step complete')
        
        if DEBUG__STEADYSTATE: break
        
