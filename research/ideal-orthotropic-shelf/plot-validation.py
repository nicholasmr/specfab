#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Validate that solving for steady-state fabric gives the correct result as forward modelled
"""

import importlib, sys, os, code # code.interact(local=locals())
import numpy as np

from lib.constants import *
from lib.plottools import *
from lib.model import *

from dolfin import *

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
from specfabpy.fenics.ice import IceFabric
import specfabpy.fenics.tools as sffenicstools

from matplotlib import rcParams, rc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
import cmcrameri.cm as cmc

FS = sfplt.setfont_tex(fontsize=11)

"""
Setup
"""

if len(sys.argv) != 3: 
    print('usage: %s [LROT|DDRX] L'%(sys.argv[0]))
    sys.exit(0)

from R24config import *
conf = R24config()
conf.rheology = 'Orthotropic'
conf.fabdyn   = str(sys.argv[1]) # fabric dynamics type [LROT|DDRX]
conf.L        = int(sys.argv[2]) # spectral truncation

lm, nlm_len = sf.init(conf.L)

tstep = 100  # time step to plot of forward model (should correspond to a steady state)
scale = 1e-3 # spatial scale

### Load model state 

model = FlowModel(1, **conf.kwargs_flowmodel())
state = load_state(conf.fname_state(tstep), conf.L, mesh=model.mesh.mesh)
model.set_state(*state) 
coords, triang, meshpts = model.mesh.get_coords(scale=scale)
fabric = model.momentum.fabric # shorthand used below
u_true = model.momentum.u # shorthand used below
pfJ = model.mesh.vertexvals_on_Q(fabric.pfJ)

### Solve for steady-state fabric

print('*** Solving for steady-state CPO field')

# Nonlinear solver does not converge unless regularization is made a bit stronger, 
#  hence LROT+DDRX fabrics cannot be perfectly recovered (unlike LROT-only fabrics).
if fabric.fabdyn=='DDRX': 
    fabric.nu_realspace *= 2

cpo = IceFabric(model.mesh.mesh, model.mesh.boundaries, L=fabric.L, modelplane=model.mesh.modelplane, symframe=fabric.symframe, \
                nu_multiplier=fabric.nu_multiplier, nu_realspace=fabric.nu_realspace)
cpo.initialize() # isotropic by default when no arguments passed
cpo.set_isotropic_BCs([DOM_ID__IN,]) # set isotropic on inflow/divide boundaries

if   fabric.fabdyn == 'LROT': Gamma0 = None
elif fabric.fabdyn == 'DDRX': Gamma0 = [fabric.Gamma0_Lilien23_lab(u_true, T) for T in np.linspace(-40+273.15,model.rheology.T,3)] # list of DDRX rate factors used to gradually approach solution 
S = sym(grad(u_true)) # assumes strain-rate and stress tensor are approximately coaxial
cpo.solvesteady(u_true, S, iota=+1, Gamma0=Gamma0, LROT_guess=True)  
pfJ_inferred = model.mesh.vertexvals_on_Q(cpo.pfJ) 

### Shared for plotting 

figscale = 1.3

kwargs_cb  = dict(pad=0.04, aspect=14, fraction=0.13, orientation='vertical')
kwargs_lbl = dict(fontsize=11+1.5, frameon=False)

x0, y0 = np.amin(meshpts[0]), np.amin(meshpts[1])
x1, y1 = np.amax(meshpts[0]), np.amax(meshpts[1])

def plot_boundaries(ax, legend=False):
    c1, c2, c3 = 'c', 'k', '#f0027f'
    kwargs = dict(lw=2.2, zorder=20, clip_on=False,)
    ax.plot([x0,x0], [y0,y1], '-', c=c1, label='Isotropic', **kwargs)
    ax.plot([x0,x1,x1,x0], [y0,y0,y1,y1], '--', c=c3, label='Free', **kwargs)
    if legend: ax.legend(loc=1, ncols=3, fancybox=False, frameon=False, columnspacing=1, bbox_to_anchor=(1.02, 1.28))

def setupaxis(ax, noy=False):
    xlims = [x0,x1]
    ylims = [y0,y1]
#    ax.axis('equal')
    ax.set_xlabel(r'$x$ ($\SI{}{\kilo\metre}$)', labelpad=-1)
    ax.set_xticks(np.arange(x0, x1+1, 20))
    ax.set_xticks(np.arange(x0, x1+1, 10), minor=True)
    ax.set_xlim(xlims)
    if not noy:
        ax.set_ylabel(r'$y$ ($\SI{}{\kilo\metre}$)')
        ax.set_yticks(np.arange(y0, y1+1, 5))
        ax.set_yticks(np.arange(y0, y1+1, 2.5), minor=True)
        ax.set_ylim(ylims)
    else:
        ax.tick_params('y', labelleft=False)

### Setup figure

figwidth = 6.3*figscale
fig = plt.figure(figsize=(figwidth, 0.65*3.6*figscale))
gs = gridspec.GridSpec(3,6, height_ratios=[0.95, -0.13, 0.5])
gs.update(left=0.065, right=0.99, top=0.90, bottom=0.015, wspace=0.3, hspace=0.8)
k = 3
ax1 = plt.subplot(gs[0,0:k])
ax2 = plt.subplot(gs[0,k:], sharey=ax1)


ROTATE_TO_XY = True

if ROTATE_TO_XY: geo, prj = sfplt.getprojection(rotation=-90, inclination=0)      
else:            geo, prj = sfplt.getprojection(rotation=-90, inclination=90)

axODFs = [plt.subplot(gs[2,ii], projection=prj) for ii in range(2*k)] # inset axes
for ax in axODFs: ax.set_global()

### pfJ 

lvls = np.arange(1,5.5+1e-3,0.5)
kwargs = dict(levels=lvls, extend='max', cmap=cmc.batlowW_r)

ax = ax1
setupaxis(ax, noy=False)
h = ax.tricontourf(triang, pfJ, **kwargs)
sfplt.panellabel(ax, 3, r'True')

ax = ax2
setupaxis(ax, noy=True)
h = ax.tricontourf(triang, pfJ_inferred, **kwargs)
sfplt.panellabel(ax, 3, r'Inferred')

for ax in (ax1,ax2):
    hcb = plt.colorbar(h, ax=ax, **kwargs_cb)
    hcb.set_label(r'$J$')
    hcb.set_ticks(lvls[0::2])
    plot_boundaries(ax, legend=(ax==ax1))

### MODFs

xi = [50, 80, 45] 
yi = [1, 6, 13]

for exp in ('true', 'infr'):

    ii0 = 0 if exp=='true' else 3
    
    for ii in np.arange(len(xi)):
        ax = axODFs[ii if exp=='true' else ii+3]
        coords = (xi[ii]/scale, yi[ii]/scale)
        nlm = (fabric.get_nlm(*coords) if exp=='true' else cpo.get_nlm(*coords))/np.sqrt(4*np.pi)
        if ROTATE_TO_XY: nlm = sf.rotate_nlm_xz2xy(nlm) # rotate to usual x--y map plane view for SSA models
        lvlset = (np.linspace(0.0, 0.4, 8), lambda x,p:'%.1f'%x) if fabric.fabdyn=='LROT' else (np.linspace(0.0, 0.3, 8), lambda x,p:'%.1f'%x) # vary top level for more contrast to see DDRX contours in horizontal plane (xy shear fabrics)
        sfplt.plotODF(nlm, cpo.lm, ax, lvlset=lvlset, showcb=False, nchunk=None)
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', negaxes=False, color='k', fontsize=FS+1.5)
        mrk = "%i"%(1+ii+ii0)
        set_subplot_label(ax, 2, r'$%s$'%(mrk), frameon=False, bbox=(-0.3,1.3))
        axp = ax1 if exp=='true' else ax2 # parent plot
        peff = [pe.withStroke(linewidth=1.5, foreground='k')]
        axp.text(xi[ii],yi[ii], mrk, fontsize=FS+2.0, color='w', path_effects=peff, ha='center', va='center', zorder=20, clip_on=False) 

    ### ODF plots legends

    ax = ax1
    x0,y0 = 175, -6.25
    ax.text(-5,   y0, r'\textit{Example of MODFs}', ha='left', va='center')

### Labels
    
bbox1=(-0.15, 1.30)
lbla, lblb = ('a','b') if fabric.fabdyn=='LROT' else ('c','d')
sfplt.panellabel(ax1, 2, r'\Large\textbf{%s}'%(lbla), bbox=bbox1, **kwargs_lbl)
sfplt.panellabel(ax2, 2, r'\Large\textbf{%s}'%(lblb), bbox=bbox1, **kwargs_lbl)

plt.savefig('validation-%s.pdf'%(fabric.fabdyn), transparent=True, dpi=200)

