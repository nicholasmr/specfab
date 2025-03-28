#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022-2024

"""
Verify FEniCS implementation of the SDM for polycrystalline olivine
by comparing with the discrete directors and FSE method

DDM = Discrete Directors Method
SDM = Spectral (continuous) Directors Method
FSE = Finite Strain Ellipsoid
"""

import copy, sys, time, code # code.interact(local=locals())
import numpy as np

from specfabpy import common as sfcom
from specfabpy import constants as sfconst
from specfabpy import discrete as sfdsc
from specfabpy.fenics.olivine import OlivineFabric
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=13)
FSLEG = FS-2
FSLBL = FS+2

from dolfin import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.tri as tri
import cartopy.crs as ccrs

#----------------
# Constants and parameters
#----------------

Nt = 60        # number of time steps to take
Ngrains = 1000 # number of DDM grains

L = 10 # spectral truncation
nu_realspace = 1e-3 # regularization magnitude (for stabilization)
nu_multiplier = 1.0

#----------------
# Initialize components
#----------------

### Mesh 

DOM_ID__BOTTOM   = 0
DOM_ID__TOP      = 1
DOM_ID__LEFT     = 2
DOM_ID__RIGHT    = 3
DOM_ID__INTERIOR = 4

def unit_mesh(XDIV=40, ZDIV=40, W=1, H=1, x0=0, z0=0):

    #----------------
    # Boundaries
    #----------------

    class BottomBoundary(SubDomain):
        def inside(self, x, on_boundary): return on_boundary and near(x[1], z0)

    class TopBoundary(SubDomain):
        def inside(self, x, on_boundary): return on_boundary and near(x[1], z0+H)

    class LeftBoundary(SubDomain):
        def inside(self, x, on_boundary): return on_boundary and near(x[0], x0)

    class RightBoundary(SubDomain):
        def inside(self, x, on_boundary): return on_boundary and near(x[0], x0+W)

    #----------------
    # Mesh
    #----------------

    mesh = RectangleMesh(Point(x0, z0), Point(x0+W, z0+H), XDIV, ZDIV) # XDIV, ZDIV, "crossed") # option "crossed" stands for crossed diagonals (number of elements=XDIV*ZDIV*4)
    cell = triangle
    norm = FacetNormal(mesh) # definition of an outer normal

    #----------------
    # Boundary partitioning
    #----------------

    boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

    boundary_parts.set_all(DOM_ID__INTERIOR) # interior of the domain

    Gamma_b = BottomBoundary()
    Gamma_b.mark(boundary_parts, DOM_ID__BOTTOM)

    Gamma_t = TopBoundary()
    Gamma_t.mark(boundary_parts, DOM_ID__TOP)

    Gamma_l = LeftBoundary()
    Gamma_l.mark(boundary_parts, DOM_ID__LEFT)

    Gamma_r = RightBoundary()
    Gamma_r.mark(boundary_parts, DOM_ID__RIGHT)

    ds = Measure("ds")[boundary_parts] 
    
    return (mesh, boundary_parts, ds, norm)

xy0 = (0,0)
xy1 = (1,1)
RES = 25

mesh, boundary_parts, ds, norm = unit_mesh(XDIV=RES, ZDIV=RES, W=xy1[0], H=xy1[1])

### Steady velocity field

Vele = VectorElement("Lagrange", mesh.ufl_cell(), 2)
V = FunctionSpace(mesh, Vele) 
v = project(Expression(("v0 + (v1-v0)/H*(pow(x[1],1))","0"), H=1, v1=1, v0=0, degree=2), V)
W = TensorFunctionSpace(mesh, 'Lagrange', 2) # for FSE field
vgrad = project(grad(v),W)

S = vgrad

### Fabric model

modelplane = 'xz'
fabric = OlivineFabric(mesh, boundary_parts, L=L, nu_realspace=nu_realspace, nu_multiplier=nu_multiplier, modelplane=modelplane, enable_SDM=True, enable_FSE=True)

### Discrete DM and FSE for validation

p = np.zeros((Nt,2)) # parcel location
p[0,:] = [0.07, 0.62]

dfse = sfdsc.DFSE()
F = np.zeros((Nt,2,2))
F[0,:,:] = dfse.F

ddm_n = sfdsc.DDM(iota=+1, N=Ngrains)
ddm_b = sfdsc.DDM(iota=-1, N=Ngrains)
n = np.zeros((Nt,2)) # largest eigenvalue vector
b = np.zeros((Nt,2)) # largest eigenvalue vector
n[0,:] = [0,0] # no preferred direction if init distribution is random
b[0,:] = [0,0] 

#----------------
# Time evolution
#----------------

t = 0.0
nn = 0

# Determine time step using CFL criterion
h_min = mesh.hmin()
v_max = abs(v.vector()[:]).max()
C_CFL = 0.5 # Courant number
dt = Constant(C_CFL*h_min/v_max)

while nn < Nt:

    nn += 1
    dt_np = float(dt)
    t += dt_np # advance time
    
    info("*** Step %i :: dt=%.2e, t=%.2e" % (nn, dt, t))

    info("-> Solving fabric evolution: continuous FEM problem")

    fabric.evolve(v, S, dt)

    info("-> Solving fabric evolution: discrete Lagrangian parcel problem")

    # Update parcel position    
    p0 = p[nn-1,:]
    vx, vy = v(*p0)
    p[nn,:] = [p0[0]+dt_np*vx, p0[1]+dt_np*vy]
    
    # Update parcel CPO
    L2 = np.reshape(vgrad(*p0), (2,2))
    L3 = sfcom.mat3d(L2, modelplane)
    ddm_b.evolve(L3, dt_np)
    ddm_n.evolve(L3, dt_np)
    dfse.evolve(L2, dt_np)
    
    # Save parcel CPO orientation
    F[nn,:,:] = dfse.F
    mi_n, eigvals_n = ddm_n.eigenframe(modelplane=modelplane)
    mi_b, eigvals_b = ddm_b.eigenframe(modelplane=modelplane)
    n[nn,:] = mi_n[0,[0,2]] # x,z components of largest eigenvalue direction
    b[nn,:] = mi_b[0,[0,2]] # x,z components of largest eigenvalue direction
        
    #---------------------------------------------
    #---------------------------------------------
        
    if nn==1 or nn%7==0:
    
        fname = 'diagnostics-convection/diagnostic-%03d.png'%(nn)
        print('[->] Plotting model state: %s'%(fname))

        cpo_n, cpo_b = fabric.SDM_n, fabric.SDM_b
        ni, bi = ddm_n.v.copy(), ddm_b.v.copy()
        pn = p[:(nn+1),:]

        discrete_interval = 7
        
        ### Figure geometry
        figscale = 0.45
        figsize = (15*figscale, 10.8*figscale)
        fig = plt.figure(figsize=figsize)
        
        ### Subplot axes
        gs = gridspec.GridSpec(1, 2, wspace=0.25, hspace=0.29, left=0.08, right=0.98, top=0.90, bottom=0.33)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])

        ### Mesh coordinates and triangulation
        coords = copy.deepcopy(mesh.coordinates().reshape((-1, 2)).T)
        triang = tri.Triangulation(*coords, triangles=mesh.cells())
        
        ### Function space for plotting
        Qele = FiniteElement("CG", mesh.ufl_cell(), 2) 
        Q = FunctionSpace(mesh, Qele)
        ex,ey = Constant((1,0)), Constant((0,1))
        
        ### S2 projection
        geo, prj = sfplt.getprojection(rotation=-90-20, inclination=50)   
                    
        ### AUX
        
        c_vel = '0.35'
        c_FSE = 'tab:brown'
        lw_FSE = 0.8
        arrscale = 21
        
        def set_axes(ax, xy0=xy0, xy1=xy1):
            dy_tick=1
            dy_tickminor=dy_tick/4
            ax.set_yticks(np.arange(xy0[1],xy1[1]+dy_tick*3,dy_tick))
            ax.set_yticks(np.arange(xy0[1],xy1[1]+dy_tick*3,dy_tickminor), minor=True)
            ax.set_ylabel('$z/H$')
            ax.set_ylim([xy0[1],xy1[1]])

            dx_tick=1
            dx_tickminor=dx_tick/4
            ax.set_xticks(np.arange(xy0[0],xy1[0]+1e-5,dx_tick))
            ax.set_xticks(np.arange(xy0[0],xy1[0]+1e-5,dx_tickminor), minor=True)
            ax.set_xlabel('$x/W$', labelpad=-1)
            ax.set_xlim([xy0[0],xy1[0]])    
        
        ### Grid
        
        L, W = xy1
        def xygrid(Nx, Ny, dN, dxfrac=1, dyfrac=1):
            dx, dy = dxfrac*L/(Nx+dN), dyfrac*W/(Ny+dN)
            x_1d = np.linspace(+dx, L-dx, Nx)
            y_1d = np.linspace(+dy, W-dy, Ny)
            xv, yv = np.meshgrid(x_1d, y_1d)            
            return xv, yv, xv.flatten(), yv.flatten()
        
        Nl1, Nl2 = 4, 5
        xv_l1, yv_l1, xf_l1, yf_l1 = xygrid(Nl1, Nl1, 2) 
        xv_l2, yv_l2, xf_l2, yf_l2 = xygrid(Nl2, Nl2, 2, dyfrac=0.5) 
           
           
        ### Discrete model
        
        ux = np.array([v.sub(0)(x,yf_l2[ii]) for ii,x in enumerate(xf_l2)]).reshape(xv_l2.shape)
        uy = np.array([v.sub(1)(x,yf_l2[ii]) for ii,x in enumerate(xf_l2)]).reshape(xv_l2.shape)
        QV = ax1.quiver(xv_l2, yv_l2, ux, uy, color=c_vel)
        
        p_list = np.array(pn)
        p_list_reduced = p_list[::discrete_interval,:]
        for ii,(pix, piy) in enumerate(p_list_reduced):
            jj = discrete_interval * ii
            eigvecs, eigvals = sfcom.eigenframe(sfcom.F2C(F[jj,:,:]), modelplane=modelplane)
            sfplt.plotFSE(ax1, (pix,piy), eigvecs, eigvals, lw=lw_FSE, ls='-', c=c_FSE)
            QV1 = ax1.quiver(pix,piy, +b[jj,0], +b[jj,1], scale=arrscale, color=sfplt.c_red)
            QV1 = ax1.quiver(pix,piy, -b[jj,0], -b[jj,1], scale=arrscale, color=sfplt.c_red)
            QV2 = ax1.quiver(pix,piy, +n[jj,0], +n[jj,1], scale=arrscale, color=sfplt.c_blue)
            QV2 = ax1.quiver(pix,piy, -n[jj,0], -n[jj,1], scale=arrscale, color=sfplt.c_blue)

        # Present CPO position
        CPOpos = (pix,piy) # p_list[-1,:] 
        CPOaxpos = (0.23, +0.01)
        kwargs_CPOpos = dict(markersize=10, markeredgewidth=1.1, markeredgecolor='k', markerfacecolor='w')
        points, = ax1.plot(*CPOpos, 'X', **kwargs_CPOpos)

        ax1.text(p_list[0,0]-0.04,p_list[0,1]-0.1, '$t=0$', fontsize=FS)
        if nn>1: ax1.text(CPOpos[0]-0.04,CPOpos[1]-0.1, r'$t=%.2f$'%(t), fontsize=FS)
           

        dx = -0.05
        plt.quiverkey(QV1, 0.58+dx, 1.075, 1.5, r'$\pm\vb{m}_1$', labelpos='E', coordinates='axes', labelsep=0.05)
        plt.quiverkey(QV2, 0.85+dx, 1.075, 1.5, r'$\pm\vb{m}_2$', labelpos='E', coordinates='axes', labelsep=0.05)
        plt.quiverkey(QV , 0.38+dx, 1.075, 0.65, r'$\vb{u}$',     labelpos='E', coordinates='axes', labelsep=0.05)

        sfplt.panellabel(ax1, 2, r'\textit{(a)}', fontsize=FSLBL, frameon=False, bbox=(-0.2,1.19))
        set_axes(ax1, xy0=xy0, xy1=xy1)
#        ax1.text(0.75, 0.05, r'$t=%.2f$'%(t), fontsize=FS)
        
        if 1:
            pos = [CPOpos,]
            CPO_plots = [{'ni':ni, 'bi':bi, 'loc':(x,y), 'axloc':(0.44+(x-0.5)/4,y*0.7+0.35), } for x,y in pos]

            for CPO in CPO_plots:

                W = 0.12 # ax width
                x_, y_ = CPOaxpos
                dxax = 0.1
                axin1 = plt.axes([x_-dxax/2,y_, W,W], projection=prj) 
                axin1.set_global()
                axin2 = plt.axes([x_+dxax/2,y_, W,W], projection=prj) 
                axin2.set_global()

                n_lat, n_colat, n_lon = sfdsc.cart2sph(CPO['ni'], deg=True)
                b_lat, b_colat, b_lon = sfdsc.cart2sph(CPO['bi'], deg=True)
                kwargs = dict(marker='o', s=0.6, linewidths=0.17, transform=geo, zorder=10)
                axin1.scatter(n_lon, n_lat, c=sfplt.c_blue, edgecolors=sfplt.c_blue, facecolors=sfplt.c_blue, **kwargs)  
                axin2.scatter(b_lon, b_lat, c=sfplt.c_red,  edgecolors=sfplt.c_red,  facecolors=sfplt.c_red,  **kwargs) 
                
                kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
                gl1 = axin1.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
                gl2 = axin2.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
                gl1.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
                gl2.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

                sfplt.plotcoordaxes(axin1, geo, axislabels='vuxi', color='k', fontsize=FSLEG)
                sfplt.plotcoordaxes(axin2, geo, axislabels='vuxi', color='k', fontsize=FSLEG)

                axin1.set_title(r'$n$', fontsize=FS)
                axin2.set_title(r'$b$', fontsize=FS)
                
            ax1.text(0.4, -0.24, 'CPO at', fontsize=FS)
            points, = ax1.plot(0.64, -0.22, 'X', clip_on=False, **kwargs_CPOpos)
        
            sfplt.panellabel(ax1, 3, r'Discrete DM and FSE', fontsize=FSLEG, frameon=True)
        
        ### Continuous model
       
        QV = ax2.quiver(xv_l2, yv_l2, ux, uy, color=c_vel)
        
        m1_n = np.zeros((Nl1,Nl1,2))
        m1_b = np.zeros((Nl1,Nl1,2))
        FSE_eigvals = np.zeros((Nl1,Nl1, 2))
        FSE_eigvecs = np.zeros((Nl1,Nl1, 2,2))
        for ii in np.arange(Nl1):
            for jj in np.arange(Nl1):
                x_, y_ = xv_l1[ii,jj], yv_l1[ii,jj]
                FSE_eigvecs[ii,jj,:,:], FSE_eigvals[ii,jj,:] = fabric.get_FSE(x_, y_) # vi = eigvecs[:,I]            
                eigvecs, eigvals = cpo_n.eigenframe(x_, y_)
                m1_n[ii,jj,:] = eigvecs[0,[0,2]] # x-z components of largest nlm eigenvalue
                eigvecs, eigvals = cpo_b.eigenframe(x_, y_)
                m1_b[ii,jj,:] = eigvecs[0,[0,2]] 

        QV1 = ax2.quiver(xv_l1, yv_l1, +m1_n[:,:,0], +m1_n[:,:,1], scale=arrscale, color=sfplt.c_blue)
        QV1 = ax2.quiver(xv_l1, yv_l1, -m1_n[:,:,0], -m1_n[:,:,1], scale=arrscale, color=sfplt.c_blue)
        QV2 = ax2.quiver(xv_l1, yv_l1, +m1_b[:,:,0], +m1_b[:,:,1], scale=arrscale, color=sfplt.c_red)
        QV2 = ax2.quiver(xv_l1, yv_l1, -m1_b[:,:,0], -m1_b[:,:,1], scale=arrscale, color=sfplt.c_red)
                
        set_axes(ax2, xy0=xy0, xy1=xy1)
        kwargsgrid = {'lw':0.075, 'color':'0.5', 'alpha':0.75}
        ax2.triplot(triang, **kwargsgrid)
        sfplt.panellabel(ax2, 2, r'\textit{(b)}', fontsize=FSLBL, frameon=False, bbox=(-0.2,1.19))
        ax2.text(0.75, 1.05, r'$t=%.2f$'%(t), fontsize=FS)
        
        for ii in np.arange(Nl1):
            for jj in np.arange(Nl1):
                x_, y_ = xv_l1[ii,jj], yv_l1[ii,jj]
                sfplt.plotFSE(ax2, (x_,y_), FSE_eigvecs[ii,jj,:,:], FSE_eigvals[ii,jj,:], lw=lw_FSE, c=c_FSE)

        if 1:
        
            points, = ax2.plot(*CPOpos, 'P', **kwargs_CPOpos)
        
            pos = [CPOpos,]
            CPO_plots = [{'nlm':cpo_n.get_nlm(x,y), 'blm': cpo_b.get_nlm(x,y), 'loc':(x,y), 'axloc':(0.44+(x-0.5)/4,y*0.7+0.35), } for x,y in pos]

            for CPO in CPO_plots:
            
                W = 0.12 # ax width
                x_, y_ = CPOaxpos
                x_ += 0.5
                dxax = 0.1
                axin1 = plt.axes([x_-dxax/2,y_, W,W], projection=prj) 
                axin1.set_global()
                axin2 = plt.axes([x_+dxax/2,y_, W,W], projection=prj) 
                axin2.set_global()
                
                sfplt.plotODF(CPO['nlm'], cpo_n.lm, axin1, cmap='Blues', lvlset=(np.linspace(0.1, 0.4, 6), lambda x,p:'%.1f'%x), showcb=False)
                sfplt.plotODF(CPO['blm'], cpo_b.lm, axin2, cmap='Reds',  lvlset=(np.linspace(0.1, 0.4, 6), lambda x,p:'%.1f'%x), showcb=False)
                sfplt.plotcoordaxes(axin1, geo, axislabels='vuxi', color='k', fontsize=FSLEG)
                sfplt.plotcoordaxes(axin2, geo, axislabels='vuxi', color='k', fontsize=FSLEG)

                axin1.set_title(r'$n$', fontsize=FS)
                axin2.set_title(r'$b$', fontsize=FS)

            ax2.text(0.4, -0.24, 'CPO at', fontsize=FS)
            points, = ax2.plot(0.64, -0.22, 'P', clip_on=False, **kwargs_CPOpos)
            
            sfplt.panellabel(ax2, 3, r'Continuous DM and FSE', fontsize=FSLEG, frameon=True)

        ### Save plot
        
        plt.savefig(fname, dpi=225)
        print('[OK] Done')
        
