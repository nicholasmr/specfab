# N. M. Rathmann <rathmann@nbi.ku.dk> and D. Lilien, 2022-2023

"""
Fabric state-space diagram quantifying the effect of truncation (L) and (hyper)regularization (Lilien et al., 2023)
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import pickle

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, '..')
from localheader import *

sys.path.append('../../')
import demolib as dl

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import integrator as sfint
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-1.5
FSANNO = FS-2.5

SELFNAME = sys.argv[0][:-3] # used as prefix for pickled files
os.system('mkdir -p specfab-state-trajectories')
def pfile(fname): return "specfab-state-trajectories/%s--%s.p"%(SELFNAME, fname) # full path name for pickled files

#--------------------
# Run options
#--------------------

DEBUG = 0 # For faster plotting (lower resolution)
            
#--------------------
# Setup
#--------------------

RESX = RESY = 1*100 if DEBUG else 15*100 

L = 8 if DEBUG else 20 # High-order truncation

Nt = 150

#--------------------
# Modeled correlations
#--------------------

def integrate_model(modtype, Nt, L, nu=1, regexpo=None):

    Gamma0 = Lambda = None
    iota, zeta = 1, 0
        
    if modtype == 'uc': mod, target = dict(type='ps', axis=2, T=+1, r=0), -0.99
    if modtype == 'ue': mod, target = dict(type='ps', axis=2, T=-1, r=0), 12

    lm, nlm_len = sf.init(L) 
    nlm, F, *_ = sfint.lagrangianparcel(sf, mod, target, Nt=Nt, iota=iota, zeta=zeta, Gamma0=Gamma0, Lambda=Lambda, nu=nu, regexpo=regexpo)
    return nlm/norm

nlm_uc20 = integrate_model('uc', Nt, 20) 
nlm_ue20 = integrate_model('ue', Nt, 20) 

nlm_uc6 = integrate_model('uc', Nt, 6)
nlm_ue6 = integrate_model('ue', Nt, 6) 

nualt = 3.5
nlm_uc6lap = integrate_model('uc', Nt, 6, nu=nualt, regexpo=1) # regexpo=1 is Laplacian diffusion
nlm_ue6lap = integrate_model('ue', Nt, 6, nu=nualt, regexpo=1) 

#--------------------
# Construct plot
#--------------------

### Setup figure

scale = 3.5
fig = plt.figure(figsize=(1*1.55*scale,1.5*scale))
gs = gridspec.GridSpec(1, 1, hspace=0)
plt.subplots_adjust(left=0.12, right=0.98, top=0.9, bottom=0.3)
ax = fig.add_subplot(gs[0,0])

xlims, ylims = [-1.45,2.65], [-1.8,3.75]
sc = np.diff(ylims)/np.diff(xlims)

mse = 7.0

### Plot enhancements

grain_params = sfconst.ice['viscoplastic']['linear'] 
(Ezz, Exz, x, y, isvalid) = dl.Eij_statemap_ice(grain_params, xlims, ylims, resx=RESX, resy=RESY, ignore_isvalid=True)

cmap = plt.cm.Blues.copy()
cmap.set_under('#d6604d')
CS = ax.contourf(x,y,Exz, levels=np.arange(0,6,1), cmap=cmap, extend='both')
hcb = plt.colorbar(CS, ax=ax, orientation='vertical') # , aspect=25
hcb.set_label(r'$E_{xz}$')

lvls = [0.01, 0.1, 0.5, 1, 1.5, 2]
CS = ax.contour(x,y,Ezz, levels=lvls, linewidths=0.7, linestyles='-', colors='k')
manual_locations = [(-0.2, 1.5)] + [(0.51, y_) for y_ in np.linspace(-1.4, 2.1, len(lvls[1:]))]
ax.clabel(CS, CS.levels, inline=True, fmt=r'$E_{zz}=%.2f$', fontsize=FS-1, manual=manual_locations)

### Model lines

h_uc20   = plot_trajectory(ax, nlm_uc20,   arrpos=9,    c=dl.c_unidir,  ls='--', endmarker=True, mse=mse)
h_ue20   = plot_trajectory(ax, nlm_ue20,   arrpos=17,   c=dl.c_planar,  ls='--', endmarker=True, mse=mse)
h_uc6    = plot_trajectory(ax, nlm_uc6,    arrpos=None, c=dl.c_unidir,  ls='-',  endmarker=True, mse=mse)
h_ue6    = plot_trajectory(ax, nlm_ue6,    arrpos=None, c=dl.c_planar,  ls='-',  endmarker=True, mse=mse)
h_uc6lap = plot_trajectory(ax, nlm_uc6lap, arrpos=None, c=dl.cl_unidir, ls=':',  endmarker=True, mse=mse)
h_ue6lap = plot_trajectory(ax, nlm_ue6lap, arrpos=None, c=dl.cl_planar, ls=':',  endmarker=True, mse=mse)

### Legend 

h_modellines = [h_ue20, h_ue6, h_ue6lap,   h_uc20, h_uc6, h_uc6lap]
legend_strings = ['Unconf. compr., L=20', 'Unconf. compr., L=6', 'Unconf. compr., L=6, no hyper dif.', \
                  'Unconf. exten., L=20', 'Unconf. exten., L=6', 'Unconf. exten., L=6, no hyper dif.']
legkwargs = {'handlelength':1.4, 'framealpha':0.7, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}
legend_modellines = ax.legend(h_modellines, legend_strings, title=r'{\bf Modelled fabric state trajectories}', title_fontsize=FSLEG, \
    loc=3, bbox_to_anchor=(-0.17,-0.5), ncol=2, fontsize=FS-1, frameon=True, **legkwargs)

#### End-member cases

dl.plot_nlm_cases(ax, FSANNO, ms=mse, show_circle=False, isolbl_above=False)

### Power spectrum cap

(n20_unidir, n40_unidir), _, _ = sfdsc.nlm_ideal_cases(norm=norm)

if 1:
    badspec = np.ones((RESY, RESX)) 
    I = np.argmin(np.abs(x-n20_unidir))
    badspec[:,I:] = 0
    I = np.argmin(np.abs(y-n40_unidir))
    badspec[I:,:] = 0
    lvls = [0, 0.5]
    plt.rcParams['hatch.color'] = sfplt.c_dgray
    plt.rcParams['hatch.linewidth'] = 0.6
    cs = ax.contourf(x, y, badspec, hatches=["\\",None], levels=lvls, cmap='gray', alpha=0)
    cs = ax.contour( x, y, badspec, levels=[0.5,], linewidths=plt.rcParams['hatch.linewidth'], colors=plt.rcParams['hatch.color'])
    ax.text(-1.2, 3.3, r'\bf Restricted', fontsize=FS-2, color='0.1')
else:
    lw = 1.1
    ax.plot([n20_unidir,n20_unidir],[-10,n40_unidir],':k', lw=lw)
    ax.plot([-10,n20_unidir],[n40_unidir,n40_unidir],':k', lw=lw)

### Misc 

ax.set_xlabel(r'$\hat{n}_2^0$')
ax.set_ylabel(r'$\hat{n}_4^0$')

ax.set_xlim(xlims)
ax.set_ylim(ylims)

# Second x axis
secax = ax.secondary_xaxis('top', functions=(n20_to_azz, azz_to_n20))
secax.set_xlabel(r'$a^{(2)}_{zz}$')

### Save figure

plt.savefig('%s.png'%(SELFNAME), dpi=175)

