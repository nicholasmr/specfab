#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

from .constants import *
from .defaultconfig import *

from specfabpy import plotting as sfplt
from specfabpy import constants as sfconst

import copy, os, sys, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp
from cmcrameri import cm as cmc

import cartopy.crs as ccrs

import cmasher as cmr
import warnings
warnings.filterwarnings("ignore")
from dolfin import *
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri
from matplotlib import rcParams, rc, colors
import matplotlib.ticker as mticker
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  

### Matplotlib font
FS = 13
rc('font',**{'family':'serif','sans-serif':['Times'],'size':FS})

### Legend style
legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'handlelength':1.9, 'labelspacing':0.25, 'columnspacing':1.25}
leglw = 0.7
gridlw = 0.5

### Default ODF plotting params
lvls_default      = np.linspace(0.0,0.4,9)
tickintvl_default = 4

def plot_diagnostics(model, conf, fname):

    print('[->] Plotting diagnostics: %s'%(fname))

    WIDTH  = model.mesh.W
    LENGTH = model.mesh.L

    #-------------------------
    # Initialize
    #-------------------------
    
    ### Figure geometry
    figscale = 0.6
#    figsize = (18*figscale, 17*figscale)
    figsize = (18*figscale, 20*figscale)
    fig = plt.figure(figsize=figsize)

    ### Subplot axes
    gs = gridspec.GridSpec(4, 4, height_ratios=[1]*4, right=1.05, wspace=0.2, hspace=-0.85)

    ax_ux  = fig.add_subplot(gs[0,:])
    ax_eig = fig.add_subplot(gs[1,:])
    ax_Eij = fig.add_subplot(gs[2,:])

    N_odfs = 4
    gsodf = gs[3,:].subgridspec(1, N_odfs, hspace=0.15, wspace=0.15)
    geo, prj = sfplt.getprojection(rotation=-90, inclination=70)   
    axODFs = [fig.add_subplot(gsodf[0,ii], projection=prj) for ii in np.arange(N_odfs)]
    for ax in axODFs: ax.set_global()
    
    ### Mesh coordiantes and triangulation
    coords, triang = model.mesh.coords, model.mesh.triang
    
    ### Function space for plotting
    Qele = FiniteElement("CG", model.mesh.mesh.ufl_cell(), 2) 
    Q = FunctionSpace(model.mesh.mesh, Qele)
    ex,ey = Constant((1,0)), Constant((0,1))
    
    #-------------------------
    # Generate plots
    #-------------------------
        
    ### VELOCITY MAP
    
    ax = ax_ux
    
    UX = project(model.momentum.u.sub(0),Q).compute_vertex_values(model.mesh.mesh)
    UY = project(model.momentum.u.sub(1),Q).compute_vertex_values(model.mesh.mesh)
    scale = Constant(MS_TO_MYR)
    umag = project(scale*sqrt(inner(model.momentum.u, model.momentum.u)),Q)
        
    # Plot contours
    ticks = 10
    plot_panel(ax, umag, model.mesh.mesh, triang, ticks, cmc.bilbao_r, r'$\vert u \vert$ (m/yr)', WIDTH,LENGTH, extend='max', dtickcbar=1)

    # Stream lines
    N = 100 
    x_1d = np.linspace(0,LENGTH,N)
    y_1d = np.linspace(0,WIDTH,N)
    xv, yv = np.meshgrid(x_1d, y_1d)
    U = griddata(coords.T, UX, (xv, yv), method='nearest')
    V = griddata(coords.T, UY, (xv, yv), method='nearest')
    start_points = np.array([(0,yi) for yi in np.arange(0.1*WIDTH,0.9*WIDTH+1e-3,5e3)])
#    ax.streamplot(xv,yv,U,V, color='0.05', linewidth=0.6, broken_streamlines=False, start_points=start_points )
    ax.streamplot(xv,yv,U,V, color='0.05', density=0.6, linewidth=0.8) 

    set_subplot_label(ax, 2, r'$\mathbf{a}$')
    dt = model.timeax[model.tt]-model.timeax[model.tt-1]
    ax.set_title(r'%s :: $n=%i$ :: $L=%i$ :: $t=%.1f yr $'%(model.rheology.rheology, model.rheology.n, model.momentum.fabric.L, model.time/MS_TO_MYR))
    
    ### Eigenvalue map
    
    ax = ax_eig
    
    deltaeig = project(abs(model.momentum.fabric.lam1-model.momentum.fabric.lam2),Q)
    lbl=r'$\vert \lambda_1 - \lambda_2 \vert$'
    ticks = np.arange(0.1, 0.7 +1e-3, 0.1)
    plot_panel(ax, deltaeig, model.mesh.mesh, triang, ticks, cmc.lapaz_r, lbl, WIDTH,LENGTH, extend='both', dtickcbar=1)

    ax.triplot(triang, lw=0.075, color='0.35')
    set_subplot_label(ax, 2, r'$\mathbf{b}$')
    
    if 1:
        lam3 = project(model.momentum.fabric.lam3, Q)
        F = lam3.compute_vertex_values(model.mesh.mesh) # vertical eigenvalue
        CS1 = ax.tricontour(triang, F, levels=[0.2, 0.4, 0.6], colors='k')
        ax.clabel(CS1, CS1.levels, inline=True, fmt='%.1f', fontsize=FS)
        h1,_ = CS1.legend_elements()
        ax.legend([h1[0],], [r'$\lambda_3$',], loc=1, framealpha=1, **legkwargs)
   
    if 1:
        Nx, Ny = 15, 7
        dx, dy = LENGTH/(Nx+1), WIDTH/(Ny+1)
        xv    = np.linspace(dx,LENGTH-dx,Nx)
        yv    = np.linspace(dy,WIDTH-dy,Ny)
        X,Y   = np.meshgrid(xv, yv)
        Xv,Yv = X.flatten(), Y.flatten()
        
        mi, _ = model.momentum.fabric.eigenframe(Xv, Yv, reduce2D=True)
        m1, m2 = mi[:,0,:], mi[:,1,:] # mi[node,i,xyz]
                
        kwargs = dict(scale=7.5**2, width=0.0025)
        QV1 = ax.quiver(Xv,Yv, +m1[:,0], +m1[:,1], color=sfplt.c_dred, **kwargs)
        QV1 = ax.quiver(Xv,Yv, -m1[:,0], -m1[:,1], color=sfplt.c_dred, **kwargs)
        QV2 = ax.quiver(Xv,Yv, +m2[:,0], +m2[:,1], color=sfplt.c_dgreen, **kwargs)
        QV2 = ax.quiver(Xv,Yv, -m2[:,0], -m2[:,1], color=sfplt.c_dgreen, **kwargs)

        plt.quiverkey(QV1, 0.85, 1.08, 1.5, r'$\pm\mathbf{m}_1$', labelpos='E', coordinates='axes', labelsep=0.05)
        plt.quiverkey(QV2, 0.95, 1.08, 1.5, r'$\pm\mathbf{m}_2$', labelpos='E', coordinates='axes', labelsep=0.05)
       
    ### Cartesian enhancements

    ax = ax_Eij
    
    F = project(model.momentum.fabric.Exx, Q)
    lbl='$E_{xx}$'
    ticks = np.arange(0.3, 1.2 +1e-3, 0.1)
    divnorm = mpl.colors.TwoSlopeNorm(vmin=ticks[0], vcenter=1, vmax=ticks[-1])
    cmap = 'BrBG'
    plot_panel(ax, F, model.mesh.mesh, triang, ticks, cmap, lbl, WIDTH,LENGTH, extend='both', dtickcbar=1, norm=divnorm)
    ax.triplot(triang, lw=0.075, color='0.35')
    
    Eij = project(model.momentum.fabric.Exz,Q).compute_vertex_values(model.mesh.mesh) # xz because that is the model horizontal plane
    CS1 = ax.tricontour(triang, Eij, levels=[1.4, 1.8], colors='k')
    ax.clabel(CS1, CS1.levels, inline=True, fmt=r'%.1f', fontsize=FS-1) # fmt=r'$E_{xx}=%.1f$'
    h1,_ = CS1.legend_elements()
    ax.legend([h1[0],], [r'$E_{xy}$',], loc=1, framealpha=1, **legkwargs)

    set_subplot_label(ax, 2, r'$\mathbf{c}$')
    
    ##############
    
    ax = ax_ux
    
    F = project(model.mesh.h, Q)
    lbl='$h$'
    ticks = np.linspace(0,model.mesh.H0,11)
    try: CS, cls = plot_panel_lines(ax_ux, F, model.mesh.mesh, triang, ticks, WIDTH,LENGTH, fs=FS+1, fmt='$%i$', rightside_up=False)
    except: pass
            

    ### ODFs
    
    if 1:
        mrk = ["%i"%(1+ii) for ii in np.arange(len(conf.odf_x))]
        for ii, ax in enumerate(axODFs):

            nlm = model.momentum.fabric.get_nlm(conf.odf_x[ii], conf.odf_y[ii])/np.sqrt(4*np.pi) 
            sfplt.plotODF(nlm, model.momentum.fabric.lm, ax, lvlset=(np.linspace(0.0, 0.4, 8), lambda x,p:'%.1f'%x), showcb=False)
            sfplt.plotcoordaxes(ax, geo, axislabels='xi', color='k', fontsize=FS+1)
            
            mi, _ = model.momentum.fabric.eigenframe(conf.odf_x[ii], conf.odf_y[ii])
            sfplt.plotmi(ax, mi, geo, marker='.', ms=10, markeredgewidth=1.0, zorder=10, colors=(sfplt.c_dred,sfplt.c_dgreen,sfplt.c_dblue))
            set_subplot_label(ax, 2, r'$({\bf %s})$'%(mrk[ii]), frameon=False, bbox=(0.0,1.2))
            
            args = (conf.odf_x[ii],conf.odf_y[ii])
            kwargs = dict(marker=r"$\mathsf{%s}$"%mrk[ii], markerfacecolor='k', markeredgecolor='k', ms=8, clip_on=False)
            kwargs2 = dict(markerfacecolor='none', markeredgecolor='k', marker='o', ms=14.5, clip_on=False)
            for ax_ in (ax_ux, ax_eig, ax_Eij):
                ax_.plot(*args, **kwargs)
                ax_.plot(*args, **kwargs2)

           
    ### SAVE PLOT
    
    gs.tight_layout(fig)
    plt.savefig(fname, dpi=130)
    
    print('[OK] Done')
    

def plot_panel(ax, z, mesh, triang, ticks, cmap, title, W1,L, dtickcbar=1, contourlvls=None, extend='max', fmt='%1.1f', cbaspect=15, ylbl='z', norm=None, locator=None):

    Z = z.compute_vertex_values(mesh)
    h = ax.tricontourf(triang, Z, levels=ticks, cmap=cmap, extend=extend, norm=norm, locator=locator)
    cbar = plt.colorbar(h, aspect=cbaspect, fraction=0.05, ax=ax)
    if locator is not None: 
        cbar.locator = locator
        cbar.minorticks_on()
    cbar.ax.set_ylabel(title)
    if dtickcbar>1:
        for label in cbar.ax.yaxis.get_ticklabels()[1::dtickcbar]: label.set_visible(False)
    setup_axes(ax, W1, L)


def plot_panel_lines(ax, z, mesh, triang, ticks, W1,L, fmt='%1.1f', lw=1.5, color='#2b8cbe', fs=FS, rightside_up=True, inline_spacing=50):

    Z = z.compute_vertex_values(mesh)
    CS = ax.tricontour(triang, Z, levels=ticks, linewidths=lw, colors=[color,])
    cls = ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=fs, rightside_up=rightside_up, inline_spacing=inline_spacing)
    setup_axes(ax, W1, L)
    return CS, cls 
        
def setup_axes(ax, W1, L):
    dy_tick=5e3
    dy_tickminor=dy_tick/2
    ymin = 0
    ymax = W1
    ax.set_yticks(np.arange(ymin,ymax+dy_tick*3,dy_tick))
    ax.set_yticks(np.arange(ymin,ymax+dy_tick*3,dy_tickminor), minor=True)
    ax.set_ylabel('$y$ (m)')
    ax.set_ylim([ymin,ymax])

    dx_tick=10e3
    dx_tickminor=dx_tick/2
    xmax = L
    ax.set_xticks(np.arange(0,xmax+1e-5,dx_tick))
    ax.set_xticks(np.arange(0,xmax+1e-5,dx_tickminor), minor=True)
    ax.set_xlabel('$x$ (m)')
   
def set_subplot_label(ax, loc, txt, frameon=True, alpha=1.0, fontsize=FS, pad=0.3, ma='none', bbox=None):
    at = AnchoredText(txt, loc=loc, prop=dict(size=fontsize), frameon=frameon, pad=pad, bbox_to_anchor=bbox, bbox_transform=ax.transAxes)
    at.patch.set_linewidth(0.7)
    ax.add_artist(at)
    
