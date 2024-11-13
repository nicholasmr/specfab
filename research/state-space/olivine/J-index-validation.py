# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2024

import copy, sys, code # code.interact(local=locals())
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

sys.path.insert(0, '..')
from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import integrator as sfint
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-1

#--------------------
# Run options
#--------------------

L = 20
Nt = 150

#--------------------
# Model integratation
#--------------------

### specfab runs

lm, nlm_len = sf.init(L) 

mod_uc, target_uc = dict(type='ps', axis=2, T=+1, r=0), -0.95
mod_ue, target_ue = dict(type='ps', axis=2, T=-1, r=0), 5

kwargs   = dict(Nt=Nt, zeta=0, nu=1, Gamma0=None, Lambda=None)
kwargs_D = dict(Nt=Nt, zeta=0, nu=1, Gamma0=None, Lambda=2.5e-2)

iota = -1 # model slip plane directions
nlm_uc1,  sf_F_uc1,  *_ = sfint.lagrangianparcel(sf, mod_uc, target_uc, iota=iota, **kwargs) 
nlm_ue1,  sf_F_ue1,  *_ = sfint.lagrangianparcel(sf, mod_ue, target_ue, iota=iota, **kwargs) 
nlmd_uc1, sfd_F_uc1, *_ = sfint.lagrangianparcel(sf, mod_uc, target_uc, iota=iota, **kwargs_D)
nlmd_ue1, sfd_F_ue1, *_ = sfint.lagrangianparcel(sf, mod_ue, target_ue, iota=iota, **kwargs_D) 

iota = +1 # model slip plane normals
nlm_uc2,  sf_F_uc2,  *_ = sfint.lagrangianparcel(sf, mod_uc, target_uc, iota=iota, **kwargs) 
nlm_ue2,  sf_F_ue2,  *_ = sfint.lagrangianparcel(sf, mod_ue, target_ue, iota=iota, **kwargs) 
nlmd_uc2, sfd_F_uc2, *_ = sfint.lagrangianparcel(sf, mod_uc, target_uc, iota=iota, **kwargs_D)
nlmd_ue2, sfd_F_ue2, *_ = sfint.lagrangianparcel(sf, mod_ue, target_ue, iota=iota, **kwargs_D) 

### Load D-Rex runs
drex_nlm_ue1, drex_F_ue1 = load_drex_run('axisymmetricExtension',   1)
drex_nlm_ue2, drex_F_ue2 = load_drex_run('axisymmetricExtension',   2)
drex_nlm_uc1, drex_F_uc1 = load_drex_run('axisymmetricCompression', 1)
drex_nlm_uc2, drex_F_uc2 = load_drex_run('axisymmetricCompression', 2)

#--------------------
# Plot
#--------------------

### Setup plot

size = 1.6
fig = plt.figure(figsize=(4*size,1.3*size))
ax = plt.subplot()
plt.subplots_adjust(left=0.08, right=1-0.01, top=0.95, bottom=0.26)

xlims = [-1, 2]
ylims = [0.75, 5]

dx, y0 = 0.035, ylims[-1]*0.42 # *0.55
c='0.0'
x0 = 0
ax.plot([x0,x0], [0,10], ':', lw=1.1, c=c)
ax.text(1*(x0-dx), y0, r'$\uparrow$ {Compression}', color=c, rotation=90, horizontalalignment='right', verticalalignment='bottom', fontsize=FSLEG)
ax.text(1*(x0+dx), y0, r'$\downarrow$ {Extension}', color=c, rotation=90, horizontalalignment='left',  verticalalignment='bottom', fontsize=FSLEG)


### Model lines

lw = 1.75
color_b = sfplt.c_red
color_n = sfplt.c_blue 

ax.plot(f_strainzz(sf_F_uc1), f_J(nlm_uc1), '-', c=color_b, lw=lw, label=r'SDM') 
ax.plot(f_strainzz(sf_F_ue1), f_J(nlm_ue1), '-', c=color_b, lw=lw) 
ax.plot(f_strainzz(sf_F_uc2), f_J(nlm_uc2), '-', c=color_n, lw=lw)
ax.plot(f_strainzz(sf_F_ue2), f_J(nlm_ue2), '-', c=color_n, lw=lw) 

ax.plot(f_strainzz(sfd_F_uc1), f_J(nlmd_uc1), '--', c=color_b, lw=lw, label=r'SDM + $\mathcal{D}$') 
ax.plot(f_strainzz(sfd_F_ue1), f_J(nlmd_ue1), '--', c=color_b, lw=lw) 
ax.plot(f_strainzz(sfd_F_uc2), f_J(nlmd_uc2), '--', c=color_n, lw=lw)
ax.plot(f_strainzz(sfd_F_ue2), f_J(nlmd_ue2), '--', c=color_n, lw=lw) 

ax.plot(f_strainzz(drex_F_uc1), f_J(drex_nlm_uc1), ':', color=color_b, lw=lw, label=r'D-Rex')
ax.plot(f_strainzz(drex_F_ue1), f_J(drex_nlm_ue1), ':', color=color_b, lw=lw)
ax.plot(f_strainzz(drex_F_uc2), f_J(drex_nlm_uc2), ':', color=color_n, lw=lw)
ax.plot(f_strainzz(drex_F_ue2), f_J(drex_nlm_ue2), ':', color=color_n, lw=lw)

### Experimental data points

#experiments = (expr_Yabe2020, expr_Boneh2014, expr_Miyazaki2013_uc, expr_Miyazaki2013_ue)
#experiments = (expr_Yabe2020,expr_Miyazaki2013_uc,expr_Miyazaki2013_ue,expr_Boneh2014)
#experiments = (expr_Boneh2014,)    
experiments = (expr_Boneh2014, expr_Yabe2020, expr_Miyazaki2013_uc, expr_Miyazaki2013_ue)
annotatelist = [1,0,0,0]

#mrk = ['s','o','o','^']
mrk = ['d','s','o','o']
#mrk = ['d',]
color_dn = '#2166ac'
color_db = '#b2182b'
color_ln = '#92c5de'
color_lb = '#f4a582'

ms = 8.0
mew = 1.4 # marker edge width

correlations = []
for expr in experiments:
    fcorr = "observed-states/%s.p"%(expr['path']) 
    print('=== Loading correlations: %s ==='%(fcorr))
    corr = pickle.load(open(fcorr, "rb"))
    correlations.append(corr) 

for axname in ['b', 'n']:
    for jj,expr in enumerate(experiments):
    
        x = correlations[jj]['%s20'%(axname)]
        y = correlations[jj]['%s40'%(axname)]
        z = correlations[jj]['%s60'%(axname)]
        J = [f_J(np.array([norm, x[ii],y[ii],z[ii]]), Iend=-1) for ii in np.arange(len(x))]

        X = np.array(expr['Fzz'])
        
        strainzz = X-1
        lbl = expr['plotname'] if axname=='n' else None #'%s ($%s$)'%(expr['plotname'], axname)
        if expr['path'] == 'Miyazaki_etal_2013__UE': lbl = None # don't show Miyazaki_etal_2013 twice in legend
        c = color_dn if (axname=='n') else color_db
        fs = 'full' if annotatelist[jj] else 'none'
        h_, = ax.plot(strainzz, J, ms=ms, markeredgewidth=mew, ls='none', color=c, fillstyle=fs, marker=mrk[jj], label=lbl, zorder=1+(len(experiments)-jj))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'handletextpad':0.6, 'borderpad':0.37, 'edgecolor':'k'}
leg = ax.legend(fontsize=FSLEG, ncol=1, frameon=False, loc='center left', bbox_to_anchor=(1, 0.5), **legkwargs)
 
[lgd.set_color('black') for lgd in leg.legendHandles]
        
legkwargs = {'handlelength':1.4, 'framealpha':0, 'fancybox':False, 'handletextpad':0.6, 'edgecolor':'k', 'fontsize':FS+1}
h1, = ax.plot([-3,-2], [0,1], '-', color=color_b)
h2, = ax.plot([-3,-2], [0,1], '-', color=color_n) 
hleg = ax.legend([h1,h2], ['$b$', '$n$'], ncol=1, loc=2, bbox_to_anchor=(0.4,1.03), labelspacing=0.15, **legkwargs)

ax.add_artist(leg)
        
### Ticks, limits, etc.

ax.set_xticks(np.arange(-1,15,1))
ax.set_xticks(np.arange(-1,15,0.5), minor=True)
ax.set_yticks(np.arange(1,20,1))
ax.set_yticks(np.arange(1,20,1), minor=True)
ax.set_xlim(xlims)
ax.set_ylim(ylims)

ax.set_ylabel(r'$J$')
ax.set_xlabel(r'$\epsilon_{zz}$')
ax.xaxis.set_label_coords(.58, -0.17)

### Parcel deformation examples

F_uc = sf_F_uc1[np.argmin(np.abs(sf_F_uc1[:,2,2] - 0.5)),:,:]
F_ue = sf_F_ue1[np.argmin(np.abs(sf_F_ue1[:,2,2] - 1.5)),:,:]

W = 0.24
ax1 = plt.axes([0.05, 0.02, W,W], projection='3d')
ax1.patch.set_visible(False)
ax2 = plt.axes([0.25, -0.01, W,W], projection='3d')
ax2.patch.set_visible(False)

kwargs = dict(azim=35, axscale=1, axislabels=True, drawinit=True, fonttex=True, lw=0.75, colorinit=sfplt.c_dred)
sfplt.plotparcel(ax1, F_uc, **kwargs)
sfplt.plotparcel(ax2, F_ue, **kwargs)

### Save fig

f = 'J-index-validation'
plt.savefig('%s.pdf'%(f), dpi=300)
#os.system('pdfcrop %s.pdf %s.pdf'%(f,f))

