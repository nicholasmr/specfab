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
from specfabpy import common as sfcom
from specfabpy.fenics.ice import IceFabric
import specfabpy.fenics.tools as sffenicstools

from matplotlib import rcParams, rc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
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
u_df = model.momentum.u # shorthand used below

# --- DEBUG ----

if 0:
    u_df = project(Expression(('7.4*x[1]','0'), degree=1), model.momentum.U)
#    u_df = project(Expression(('3.4*x[0]','0'), degree=1), model.momentum.U)
    
    D2 = sym(grad(u_df))
    D = sfcom.mat3d(D2, fabric.modelplane)

    x,y,z = np.eye(3)
    q = z
#    q = (z+x)
    q /= np.linalg.norm(q)
    rnlm = sf.nlm_to_rnlm(sf.nlm_ideal(q, 0, conf.L), fabric.nlm_len)
    fabric.initialize(wr=rnlm)

### Calculate fields 

Q0 = FunctionSpace(model.mesh.mesh, "DG", 0)

shearfrac_df, *_ = model.momentum.fabric.enhancementfactor.shearfrac_SSA(u_df)
shearfrac        = shearfrac_df.compute_vertex_values(model.mesh.mesh)

chi_df = model.momentum.fabric.enhancementfactor.chi(u_df, fabric.s, verbose=True)
chi    = chi_df.compute_vertex_values(model.mesh.mesh)
print(chi[-10:])

### Shared for plotting 

figscale = 1.3

kwargs_cb  = dict(pad=0.04, aspect=14, fraction=0.13, orientation='vertical')
kwargs_lbl = dict(fontsize=11+2.5, frameon=False)

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

fig = plt.figure(figsize=(6*figscale, 0.4*3.6*figscale))
gs = gridspec.GridSpec(1,2, height_ratios=[1])
gs.update(left=0.065, right=0.97, top=0.85, bottom=0.2, wspace=0.15)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1], sharey=ax1)

### Plot fields

lvls1 = np.arange(0,1+1e-3,0.1)
lvls2 = np.arange(0,1+1e-3,0.1)

ax = ax1
setupaxis(ax, noy=False)
h = ax.tricontourf(triang, shearfrac, cmap='Spectral_r', extend='neither', levels=lvls1)
hcb = plt.colorbar(h, ax=ax, **kwargs_cb)
#hcb.set_label(r'{\small stretching} \quad\;$\gamma$\;\quad {\small shearing}')
hcb.set_label(r'{\fontsize{10}{10}\selectfont $\leftarrow$ stretching} \;\; $\gamma$\;\; {\fontsize{10}{10}\selectfont shearing $\rightarrow$}')
plot_boundaries(ax, legend=True)
    
ax = ax2
setupaxis(ax, noy=True)
h = ax.tricontourf(triang, chi, cmap='PiYG', extend='neither', levels=lvls2)
hcb = plt.colorbar(h, ax=ax, **kwargs_cb)
#hcb.set_label(r'{\small incompatible} \quad$\chi$\quad {\small compatible}')
hcb.set_label(r'{\fontsize{10}{10}\selectfont $\leftarrow$ incomp.} \;\; $\chi$\;\; {\fontsize{10}{10}\selectfont compatible $\rightarrow$}')
hcb.set_ticks(lvls2[0::2])
plot_boundaries(ax, legend=False)

plt.savefig('biaseval-%s.pdf'%(fabric.fabdyn), transparent=True, dpi=200)

