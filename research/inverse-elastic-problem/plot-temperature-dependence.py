# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Plot grain elastic parameters as a function of temperature, including the inferred parameter values.
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

from Cij import * # elastic constants
from plottools import *

# Use alternative labels for RSPA article?
if 1:
    Cij_exprname = {\
        'Bennett1968' : r'Bennett [39]', \
        'Gammon1983'  : r'Gammon \textit{et al.} [6]', \
        'Gagnon1988'  : r'Gagnon \textit{et al.} [65]', \
        'Dantl1968'   : r'Dantl [66]', \
        'Brockamp1964': r'Brockamp \& Querfurth [67]', \
        'Bass1957'    : r'Bass \textit{et al.} [68]', \
        'Green1956'   : r'Green \& Mackinnon [69]', \
    }

### Figure setup

scale = 3.45
fig = plt.figure(figsize=(2.0*scale,2*1.05*scale))
gs = gridspec.GridSpec(2, 3)
gs.update(hspace=0.425, wspace=0.4, left=0.085, right=0.988, top=0.955, bottom=0.350)

ax_lam  = fig.add_subplot(gs[0, 0])
ax_mu   = fig.add_subplot(gs[0, 1])
ax_gam  = fig.add_subplot(gs[0, 2])
ax_Elam = fig.add_subplot(gs[1, 0])
ax_Emu  = fig.add_subplot(gs[1, 1])
ax_Egam = fig.add_subplot(gs[1, 2])

axlist = [ax_lam, ax_mu, ax_gam,  ax_Elam, ax_Emu, ax_Egam]

### Style

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'columnspacing':0.5, 'edgecolor':'k', \
            'loc':"upper left", 'borderaxespad':0, 'fontsize':FS-1}

ms, mew = 7.5, 1.4 # markersize, markeredgewidth
xlims = [-25,-10+1] # temperature range on x axis

legbboxy = -0.37 # y displacement of bbox for legends
lwleg = 0.8 # linewidth of legend box 

color = {\
    'Bennett1968' : 'tab:green',\
    'Gammon1983'  : 'tab:orange',\
    'Gagnon1988'  : 'tab:blue',\
    'Dantl1968'   : 'tab:red',\
    'Brockamp1964': 'tab:brown',\
    'Bass1957'    : 'tab:pink',\
    'Green1956'   : 'tab:purple',\
}

### Init

allexperiments = Cij_exprname.keys()

### Temperature relations

hlist   = []
hlabels = []

def plot_temperature_relation(expr):
    T_range = np.linspace(xlims[0],xlims[1],20)
    def nomvec(gv,jj): return np.array([gv[ii,jj].n for ii in range(len(T_range))])
    def varvec(gv,jj): return np.array([gv[ii,jj].n+gv[ii,jj].s*np.array([-1,1]) for ii in range(len(T_range))])
    g = np.array([Cij_to_g(Cij_Tfunc[expr](T)) for T in T_range])
    g[:,:3] *= 1e-9
    kwargs = {'color':color[expr], 'lw':1.1, 'ls':'-'}
    kwargsfill = {'color':color[expr], 'ec':'none', 'alpha':0.25}
    for pp, ax in enumerate(axlist):
        ax.fill_between(T_range, varvec(g,pp)[:,0], varvec(g,pp)[:,1], **kwargsfill)
        p1 = ax.plot(T_range, nomvec(g,pp), **kwargs)
        p2 = ax.fill(np.NaN, np.NaN, **kwargsfill)
        if pp == 3:
            hlist.append((p2[0], p1[0]))
            hlabels.append(Cij_exprname[expr])

plot_temperature_relation('Gagnon1988')
plot_temperature_relation('Gammon1983')
plot_temperature_relation('Dantl1968')

### Single temperature experiments

for ii, expr in enumerate(allexperiments):
    if expr in ['Gammon1983','Dantl1968','Gagnon1988']: continue # skip because we are plotting the temperature relations instead
    g = Cij_to_g(Cij_ufloat[expr])
    g[0:3] *= 1e-9 # to GPa
    kwargs = {'marker':'o', 'color':color[expr], 'ls':'none', 'ms':ms, 'fillstyle':'none', 'markeredgewidth':mew, 'capsize':3}
    for pp, ax in enumerate(axlist):
        heb = ax.errorbar(Cij_T[expr], g[pp].n, g[pp].s, **kwargs)
        if pp == 3:
            hlist.append(heb)
            hlabels.append(Cij_exprname[expr])

hleg = ax_Elam.legend(hlist, hlabels, title=r'literature estimates', bbox_to_anchor=(-0.20, legbboxy), **legkwargs)
hleg.get_frame().set_linewidth(lwleg)
hleg.get_title().set_fontsize(FSSMALL)
ax_Elam.add_artist(hleg)

#### Inferred values using Thomas et al. (2021) and Lutz et al. (2022) experimental data

cred   = '#e41a1c'
cblue  = '#377eb8'
cbrown = '#a65628'

mrklist = ['^']*3 + ['<']*3 + ['v']*3
clist   = [cred, cblue, cbrown]*3 

allinferred = g_inferred.keys()
T = -23 # Temperature of Lutz et al. (2022) experiments
h_inferred = [0]*9
labels = ['']*9
for ii, expr in enumerate(allinferred):
    g_ = g_inferred[expr]
    g = [g_[0],g_[1],g_[0]+2*g_[1], g_[2],g_[3],g_[4]]
    kwargs = {'marker':mrklist[ii], 'color':clist[ii], 'ls':'none', 'fillstyle':'none', 'ms':ms, 'markeredgewidth':mew}
    for pp, ax in enumerate(axlist):
        h, = ax.plot(T, g[pp], **kwargs)
    h_inferred[ii] = h
    labels[ii] = g_exprname[expr]

hleg = ax_Emu.legend(handles=h_inferred, labels=labels, title=r'inferred from ice-core samples', bbox_to_anchor=(+0.05, legbboxy), ncol=3, **legkwargs)
hleg.get_frame().set_linewidth(lwleg)
hleg.get_title().set_fontsize(FSSMALL)

### Axis setup etc.
   
def setup_axes(ax, label, ylabel='', ylims=None, yticks=None, yticks_minor=None):
    ax.set_xlabel(r'$T$ (\SI{}{\degreeCelsius})')
    ax.set_ylabel(ylabel)
    
    ylims = ax.get_ylim()
    if yticks is not None: ax.set_yticks(yticks)
    if yticks_minor is not None: ax.set_yticks(yticks_minor, minor=True)
    ax.set_ylim(ylims)
    if ylims is not None: ax.set_ylim(ylims)
        
    ax.set_xticks(np.arange(-50,10,10))
    ax.set_xticks(np.arange(-50,10,5), minor=True)
    ax.set_xlim(xlims)
    
#    writeSubplotLabel(ax, 2, r'\bf{%s}'%label, fontsize=FS)
    writeSubplotLabel(ax, 2, r'(\textit{%s})'%label, fontsize=FS+0.5, frameon=False, alpha=0.0, pad=0.0, bbox=(-0.3,1.2))

yticks       = np.arange(0,20,1)
yticks_minor = np.arange(0,20,1/4)
setup_axes(ax_lam,  'a', yticks=yticks, yticks_minor=yticks_minor, ylabel=r'$\lambda$ (GPa)')
setup_axes(ax_mu,   'b', yticks=yticks, yticks_minor=yticks_minor, ylabel=r'$\mu$ (GPa)')
setup_axes(ax_gam,  'c', yticks=yticks, yticks_minor=yticks_minor, ylabel=r'$\gamma$ (GPa)')
yticks       = np.arange(0,2,0.1)
yticks_minor = np.arange(0,2,0.1/4)
setup_axes(ax_Elam, 'd', yticks=yticks, yticks_minor=yticks_minor, ylabel=r'$\hat{\lambda}$')
setup_axes(ax_Emu,  'e', yticks=yticks, yticks_minor=yticks_minor, ylabel=r'$\hat{\mu}$')
setup_axes(ax_Egam, 'f', yticks=yticks, yticks_minor=yticks_minor, ylabel=r'$\hat{\gamma}$')

### Save figure

fname = 'plots/temperature-dependence.png'
#fname = 'plots/temperature-dependence.pdf'
print('** Saving %s'%(fname))
plt.savefig(fname, dpi=300)

### AUX / DEBUG

if 0:

    gB68,T = Cij_to_g(Cij['Bennett1968']), Cij_T['Bennett1968'] 
    (lam,mu,gam,Elam,Emu,Egam) = gB68
#    print(gB68)

    rhoi = 917 # density of ice
    PSmap = lambda x: np.sqrt(x/rhoi) # P/S wave velocity
    
    # Analytical phase velocities for an isotropic fabric
    vS_est = PSmap( 1/15*lam*(1+1*Egam-2*Elam) + 1/15*mu*(7+2*Egam+6*Emu) )
    vP_est = PSmap( 1/15*lam*(8+3*Egam+4*Elam) + 2/15*mu*(8+3*Egam+4*Emu) )

    # Calculated phase velocities for isotropic fabric (using specfab)
    from specfabpy import specfabpy as sf 
    lm, nlm_len = sf.init(4) # L=4 suffices for this purpose
    nlm = np.zeros(nlm_len)
    nlm[0] = 1/np.sqrt(4*np.pi)
    alpha = 1 # Reuss-Voigt weight
    theta, phi = 0,0 # propagation direction (doesn't matter, considering an isotropic fabric)
    vi = sf.Vi_elastic_tranisotropic(nlm, alpha, gB68, rhoi, theta,phi)
    vP,vS1,vS2  = vi[2,0], vi[0,0], vi[1,0]

    # Compare...
    print(vS_est,vS_est, vS1,vS2)
    print(vP_est, vP)

        
