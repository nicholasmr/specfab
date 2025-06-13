#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Plot idealized model results for inferred equivalent isotropic enhancement factor
"""

import importlib, sys, os, code # code.interact(local=locals())
import numpy as np

from lib.constants import *
from lib.plottools import *
from lib.model import *

from dolfin import *

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
from specfabpy import constants as sfconst

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

if len(sys.argv) != 4: 
    print('usage: %s [LROT|DDRX] L nt'%(sys.argv[0]))
    sys.exit(0)

from R25config import *
conf = R25config()
conf.rheology = 'Orthotropic'
conf.fabdyn   = str(sys.argv[1]) # fabric dynamics type [LROT|DDRX]
conf.L        = int(sys.argv[2]) # spectral truncation
tstep         = int(sys.argv[3]) # time step to plot

fext = 'pdf'
#fext = 'png'
fsave = '%s/stateplots/%i'%(conf.path_output(), tstep)
#print(fsave)

lm, nlm_len = sf.init(conf.L)

#tstep = 100  # time step to plot of forward model (should correspond to a steady state)
scale = 1e-3 # spatial scale

Emin, Emax = sfconst.ice['viscoplastic']['CAFFEdagger'] # Emin (Ec), Emax (Es) effective values for CAFFE to match EIE

### Load model state 

model = FlowModel(1, **conf.kwargs_flowmodel())
state = load_state(conf.fname_state(tstep), conf.L, mesh=model.mesh.mesh)
model.set_state(*state) 
coords, triang, meshpts = model.mesh.get_coords(scale=scale)
h    = model.mesh.vertexvals_on_Q(model.mesh.h)
ux   = model.mesh.vertexvals_on_Q(model.momentum.u.sub(0))
uy   = model.mesh.vertexvals_on_Q(model.momentum.u.sub(1))
fabric = model.momentum.fabric
Exij = [model.mesh.vertexvals_on_Q(fabric.Exij[kk]) for kk in range(6)]
pfJ  = model.mesh.vertexvals_on_Q(fabric.pfJ)

### Shared for plotting 

figscale = 1.3

E0, E1, dE = (0.5, 1.5, 0.1) if conf.fabdyn=='LROT' else (0.5, 3.0, 0.1)
lvls = np.arange(E0, E1+1e-3, dE)
divnorm=colors.TwoSlopeNorm(vmin=np.amin(lvls), vcenter=1, vmax=np.amax(lvls))
kwargs_E = dict(levels=lvls, norm=divnorm, extend='both', cmap='PuOr_r')
E_ticks = lvls[0::5] if conf.fabdyn=='LROT' else [0.5,1,2,3]

kwargs_cb = dict(pad=0.04, aspect=14, fraction=0.13, orientation='vertical')
kwargs_lbl = dict(fontsize=11+1.5, frameon=False)

x0, y0 = np.amin(meshpts[0]), np.amin(meshpts[1])
x1, y1 = np.amax(meshpts[0]), np.amax(meshpts[1])

def plot_boundaries(ax, legend=False):
    c1, c2, c3 = 'c', 'k', '#f0027f'
    kwargs = dict(lw=2.2, zorder=20, clip_on=False,)
    ax.plot([x0,x0], [y0,y1], '-', c=c1, label='Isotropic', **kwargs)
    ax.plot([x0,x1,x1,x0], [y0,y0,y1,y1], '--', c=c3, label='Free', **kwargs)
    if legend: ax.legend(loc=1, ncols=3, fancybox=False, frameon=False, columnspacing=1, bbox_to_anchor=(1.02, 1.305))

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

"""
Orthotropic steady state
"""

### Setup figure

figwidth = 6.3*figscale
fig = plt.figure(figsize=(figwidth, 3.6*figscale))
gs = gridspec.GridSpec(4,6, height_ratios=[1, 1, -0.2, 0.6])
gs.update(left=0.065, right=0.98, top=0.94, bottom=0.015, wspace=0.3, hspace=0.8)
k = 3
ax1 = plt.subplot(gs[0,0:k])
ax2 = plt.subplot(gs[0,k:], sharey=ax1)
ax3 = plt.subplot(gs[1,0:k])
ax4 = plt.subplot(gs[1,k:], sharey=ax3)

ROTATE_TO_XY = True

if ROTATE_TO_XY: geo, prj = sfplt.getprojection(rotation=-90, inclination=0)      
else:            geo, prj = sfplt.getprojection(rotation=-90, inclination=90)

axODFs = [plt.subplot(gs[3,ii], projection=prj) for ii in range(2*k)] # inset axes
for ax in axODFs: ax.set_global()

velscale = 3.17098e-8 # m/yr to m/s

### Velocity 

ax = ax1
setupaxis(ax, noy=False)
umag = np.sqrt(np.power(ux,2) + np.power(uy,2))
umag *= 1/velscale # m/yr to m/s
lvlsu = np.arange(2,27+1,3) if conf.fabdyn == 'LROT' else np.arange(4,44+1,4)
cs = ax.tricontourf(triang, umag, levels=lvlsu, extend='both', cmap=cmc.lipari)
hcb = plt.colorbar(cs, **kwargs_cb)
hcb.set_label(r'%s (\SI{}{\metre\per\year})'%('u'))
hcb.set_ticks(lvlsu[0::2])
#ax.triplot(triang, **{'lw':0.095, 'color':'0.5', 'alpha':0.8})
    
lbl='$h$'
ticks = np.linspace(0, 500, 11)
CS = ax.tricontour(triang, h, levels=ticks, linewidths=1.5, colors=['w',])
cls = ax.clabel(CS, CS.levels, inline=True, fmt='$h=%i\,$m', rightside_up=True, inline_spacing=8)

plot_boundaries(ax, legend=True)

### pole figure J

ax = ax2
setupaxis(ax, noy=True)
lvls = np.arange(1,5.5+1e-3,0.5)
cs = ax.tricontourf(triang, pfJ, levels=lvls, extend='max', cmap=cmc.batlowW_r)
hcb = plt.colorbar(cs, **kwargs_cb)
hcb.set_label(r'$J$')
hcb.set_ticks(lvls[0::2])

plot_boundaries(ax, legend=False)

### Exx

ax = ax4
setupaxis(ax, noy=True)
h = ax.tricontourf(triang, Exij[0], **kwargs_E)
hcb = plt.colorbar(h, ax=ax, **kwargs_cb)
hcb.set_ticks(E_ticks)
hcb.set_label(r'$E_{xx}$')

plot_boundaries(ax)

### Exy

ax = ax3
setupaxis(ax, noy=False)
h = ax.tricontourf(triang, Exij[4], **kwargs_E) # model frame is xz, so this is the xy enhancement in plotted cordinate system
hcb = plt.colorbar(h, ax=ax, **kwargs_cb)
hcb.set_ticks(E_ticks)
hcb.set_label(r'$E_{xy}$')

plot_boundaries(ax)

### CPO examples

#xi = [50, 50, 80, 50, 60, 80] 
#yi = [1, 5, 6, 10, 14, 10]

xi = [50, 50, 80, 50, 60, 80] 
yi = [1, 5, 4, 10, 14, 12]

if 1:
    
    for ii in np.arange(len(axODFs)):
        ax = axODFs[ii]
        nlm = fabric.get_nlm(xi[ii]/scale, yi[ii]/scale)/np.sqrt(4*np.pi) 
        if ROTATE_TO_XY: nlm = sf.rotate_nlm_xz2xy(nlm) # rotate to usual x--y map plane view for SSA models
        lvlset = (np.linspace(0.0, 0.4, 8), lambda x,p:'%.1f'%x) if conf.fabdyn=='LROT' else (np.linspace(0.0, 0.3, 8), lambda x,p:'%.1f'%x) # vary top level for more contrast to see DDRX contours in horizontal plane (xy shear fabrics)
        sfplt.plotODF(nlm, fabric.lm, ax, lvlset=lvlset, showcb=False, nchunk=None)
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', negaxes=False, color='k', fontsize=FS+1.5)
        mi, _ = sfcom.eigenframe(nlm, symframe=fabric.symframe, modelplane='xy') 
        sfplt.plotmi(ax, mi, geo)
        
        num = "%i"%(1+ii)
        ax.add_artist(AnchoredText(num, loc=2, prop=dict(color='k', size=FS+1.5), frameon=False, pad=0.3, bbox_to_anchor=(-0.3,1.3), bbox_transform=ax.transAxes))
        peff = [pe.withStroke(linewidth=1.5, foreground='k')]
        for ax_ in (ax2,ax3,ax4): 
            ax_.text(xi[ii],yi[ii], num, fontsize=FS+2.0, color='w', path_effects=peff, ha='center', va='center', zorder=20, clip_on=False) 

### ODF plots legends

ax = ax3
ax.text(-5, -7.6, r'\textit{Example of MODFs}', ha='left', va='center')
X0, Y0 = 31.5, -11
ax.scatter([X0-3.2,],[Y0-0.5,], marker='o', s=16, c=sfplt.c_dred, clip_on=False)
ax.text(X0, Y0-0.6, r'$\vb{m}_i$', ha='left', va='center')

### Labels
    
bbox1=(-0.15, 1.36)
lblprefix = r'\fontsize{13}{13}\selectfont\textbf'
sfplt.panellabel(ax1, 2, r'%s{a}'%(lblprefix), bbox=bbox1, **kwargs_lbl)
sfplt.panellabel(ax2, 2, r'%s{b}'%(lblprefix), bbox=bbox1, **kwargs_lbl)
sfplt.panellabel(ax3, 2, r'%s{c}'%(lblprefix), bbox=bbox1, **kwargs_lbl)
sfplt.panellabel(ax4, 2, r'%s{d}'%(lblprefix), bbox=bbox1, **kwargs_lbl)
    
#plt.savefig('steadystate-%s.pdf'%(conf.fabdyn), transparent=True, dpi=200)
plt.savefig('%s-state.%s'%(fsave, fext), transparent=True, dpi=200)

"""
Isotropic enhancements
"""

conf.rheology = 'Isotropic'
isomodel = FlowModel(1, **conf.kwargs_flowmodel())
conf.rheology = 'Orthotropic' # swtich back to laod the orthotropic state into the isotropic model
state = load_state(conf.fname_state(tstep), conf.L, mesh=isomodel.mesh.mesh)
isomodel.set_state(*state) 

def plot_E(E, u, u_true, u_noE, lbla, lblb, fname, lvlsdu=[], labelpos=None):

    ### Setup figure

    fig = plt.figure(figsize=(figwidth, 1.375*figscale))
    gs = gridspec.GridSpec(1,2)
    gs.update(left=0.065, right=0.97, top=0.83, bottom=0.21, wspace=0.09, hspace=0.1)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1], sharey=ax1)

    ### Plot E
        
    ax = ax1
    setupaxis(ax, noy=False)
    h = ax.tricontourf(triang, E, **kwargs_E)
    hcb = plt.colorbar(h, ax=ax, **kwargs_cb)
    hcb.set_ticks(E_ticks)
    hcb.set_label(r'$E$')
    plot_boundaries(ax)
    
    #peff = [pe.withStroke(linewidth=1.5, foreground='k')]
    #kwargs = dict(fontsize=FS+2, color='w', path_effects=peff, ha='center', va='center',)
    kwargs = dict(fontsize=FS+2, color='k', ha='center', va='center')
    h1 = ax.text(55, 11.5, r'$E\rightarrow E_{xy}$', **kwargs)
    h2 = ax.text(55, 3, r'$E\rightarrow E_{xx}$', **kwargs)

    ### Plot velocity missfit
    
    ax = ax2
    setupaxis(ax, noy=True)
    lvlsu = np.arange(-15,15+1,2.5)
    divnorm=colors.TwoSlopeNorm(vmin=np.amin(lvlsu), vcenter=0, vmax=np.amax(lvlsu))
    kwargs_misfit = dict(levels=lvlsu, norm=divnorm, extend='both', cmap='RdBu_r')
    du = 100*np.divide(u-u_true, u_true)
    h = ax.tricontourf(triang, du, **kwargs_misfit)
    hcb = plt.colorbar(h, ax=ax, **kwargs_cb)
    hcb.set_ticks(lvlsu[0::2])
#    hcb.set_label(r'$\Delta u/u$ (\%)')
    hcb.set_label(r'Velocity error (\%)')

    du = 100*np.divide(u_noE-u_true, u_true)
    CS = ax.tricontour(triang, du, levels=lvlsdu, linestyles='--', colors='k')
    ax.clabel(CS, CS.levels, inline=True, fmt=r'\SI{%i}{\percent}', manual=labelpos, rightside_up=True, inline_spacing=6)

    nm, lbl = CS.legend_elements()
    ax.legend([nm[0],], ['Velocity error for $E=1$ (Glen)'], loc=1, fancybox=False, frameon=False, bbox_to_anchor=(1.04,1.31))
    
    plot_boundaries(ax)

    ### Save

    bbox1=(-0.15, 1.375)
    sfplt.panellabel(ax1, 2, r'\Large\textbf{%s}'%(lbla), bbox=bbox1, **kwargs_lbl)
    sfplt.panellabel(ax2, 2, r'\Large\textbf{%s}'%(lblb), bbox=bbox1, **kwargs_lbl)
    
#    plt.savefig('steadystate-%s-%s.pdf'%(conf.fabdyn,fname), transparent=True, dpi=200)
    plt.savefig('%s-%s.%s'%(fsave, fname, fext), transparent=True, dpi=200)


#    u_guess = interpolate(Constant((1e-8,0)), isomodel.momentum.U)
#    kwargs_integrate = dict(u_guess=None, tolmul=1e-12)
kwargs_integrate = dict(u_guess=None, relaxation=0.45, tolmul=5e-9)

vecmag_on_Q = lambda X: isomodel.mesh.vertexvals_on_Q(sqrt(dot(X,X)))
u_true_df = model.momentum.u
u_true = vecmag_on_Q(u_true_df)

if 1:

    ### EIE

    E_EIE_df = isomodel.momentum.set_E_EIE(u_true_df, fabric.Eij, fabric.mi, model.momentum.rheology.n) # note arguments are from the orthotropic model!
    E_EIE = isomodel.mesh.vertexvals_on_Q(E_EIE_df)
    isomodel.tt = 0 # force fine-tolerance
    isomodel.integrate(0, **kwargs_integrate)
    u_EIE = vecmag_on_Q(isomodel.momentum.u)

    _ = isomodel.momentum.set_E_EIE(u_true_df, [Constant(1)]*6, fabric.mi, model.momentum.rheology.n) # note arguments are from the orthotropic model!
    isomodel.tt = 0 # force fine-tolerance
    isomodel.integrate(0, **kwargs_integrate)
    u_noE = vecmag_on_Q(isomodel.momentum.u)
    lvlsdu   = [-30,-20,-10,0,10] if conf.fabdyn == 'LROT' else [-60, -50, -40, -30]
    labelpos = [(xi,5) for xi in (25,60,80,88,98)] if conf.fabdyn == 'LROT' else  [(xi,5) for xi in (20,60,85)]

    plot_E(E_EIE, u_EIE, u_true, u_noE, 'e', 'f', 'EIE', lvlsdu=lvlsdu, labelpos=labelpos)

    ### CAFFE 

    E_CAFFE_df = isomodel.momentum.set_E_CAFFE(fabric.s, u_true_df, Emin, Emax, 1) # note arguments are from the orthotropic model!
    E_CAFFE = isomodel.mesh.vertexvals_on_Q(E_CAFFE_df)
    isomodel.tt = 0 # force fine-tolerance
    isomodel.integrate(0, **kwargs_integrate)
    u_CAFFE = vecmag_on_Q(isomodel.momentum.u)

    plot_E(E_CAFFE, u_CAFFE, u_true, u_noE, 'g', 'h', 'CAFFE', lvlsdu=lvlsdu, labelpos=labelpos)
    
"""
Velocity error for SI
"""

def plot_velcmp(ulist, names, fname='velcmp', lvlsdu=[], labelpos=None):

    fig = plt.figure(figsize=(figwidth, 2*1.4*figscale))
    gs = gridspec.GridSpec(2,2)
    gs.update(left=0.065, right=0.98, top=0.91, bottom=0.11, wspace=0.2, hspace=0.6)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1], sharey=ax1)
    ax3 = plt.subplot(gs[1,0])
    ax4 = plt.subplot(gs[1,1], sharey=ax3)

    for ii, ax in enumerate((ax1,ax2,ax3,ax4)):
        setupaxis(ax, noy=False)
        cs = ax.tricontourf(triang, ulist[ii]/velscale, levels=lvlsu, extend='both', cmap=cmc.lipari)
        hcb = plt.colorbar(cs, **kwargs_cb)
        hcb.set_label(r'%s (\SI{}{\metre\per\year})'%('u'))
        hcb.set_ticks(lvlsu[0::2])
        plot_boundaries(ax)
        sfplt.panellabel(ax, 2, r'\textbf{%s}'%(chr(ord('a')+ii)), bbox=(-0.07, 1.35), **kwargs_lbl)
        sfplt.panellabel(ax, 9, names[ii], bbox=(0.5, 1.32), **kwargs_lbl)

    plt.savefig('%s-%s.pdf'%(fsave, fname), transparent=True, dpi=200)

if 1:
    ulist = (u_true, u_noE, u_EIE, u_CAFFE)
#    ulist = [u_true]*4 # debug
    names = ('Orthotropic rheology', 'Glen ($E=1$)', 'Glen w/ EIE', 'Glen w/ CAFFE')
    plot_velcmp(ulist, names)

