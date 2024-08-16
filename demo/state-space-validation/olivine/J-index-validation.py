# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

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
fig = plt.figure(figsize=(3.3*size,1.8*size))
ax = plt.subplot()
k=.11
plt.subplots_adjust(left=k, right=1-0.02, top=0.95, bottom=0.28)

xlims = [-1,3]
ylims = [0.75,6]

# Background patches
#rect1 = plt.Rectangle((0,ylims[0]), 1, 20, color='#c7eae5')
#rect2 = plt.Rectangle((1,ylims[0]), 20, 20, color='#f6e8c3')
#ax.add_patch(rect1)
#ax.add_patch(rect2)
#dx, y0 = 0.035, ylims[-1]*0.20 # *0.55
#ax.text(1*(1-dx), y0, r'{\bf Compression}', color='#01665e', rotation=90, horizontalalignment='right', verticalalignment='bottom', fontsize=FSLEG)
#ax.text(1*(1+dx), y0, r'{\bf Extension}',   color='#8c510a', rotation=90, horizontalalignment='left', verticalalignment='bottom', fontsize=FSLEG)

dx, y0 = 0.035, ylims[-1]*0.24 # *0.55
c='0.0'
x0 = 0
ax.plot([x0,x0], [0,10], ':', lw=1.1, c=c)
ax.text(1*(x0-dx), y0, r'$\uparrow$ {Compression}', color=c, rotation=90, horizontalalignment='right', verticalalignment='bottom', fontsize=FSLEG)
ax.text(1*(x0+dx), y0, r'$\downarrow$ {Extension}', color=c, rotation=90, horizontalalignment='left',  verticalalignment='bottom', fontsize=FSLEG)


### Model lines

lw = 1.75
color_b = sfplt.c_red
color_n = sfplt.c_blue 

ax.plot(f_strainzz(sf_F_uc1), f_J(nlm_uc1), '-', c=color_b, lw=lw, label=r'SDM ($b$)') 
ax.plot(f_strainzz(sf_F_ue1), f_J(nlm_ue1), '-', c=color_b, lw=lw) 
ax.plot(f_strainzz(sf_F_uc2), f_J(nlm_uc2), '-', c=color_n, lw=lw, label=r'SDM ($n$)') 
ax.plot(f_strainzz(sf_F_ue2), f_J(nlm_ue2), '-', c=color_n, lw=lw) 

ax.plot(f_strainzz(sfd_F_uc1), f_J(nlmd_uc1), '--', c=color_b, lw=lw, label=r'SDM + $\mathcal{D}$ ($b$)') 
ax.plot(f_strainzz(sfd_F_ue1), f_J(nlmd_ue1), '--', c=color_b, lw=lw) 
ax.plot(f_strainzz(sfd_F_uc2), f_J(nlmd_uc2), '--', c=color_n, lw=lw, label=r'SDM + $\mathcal{D}$ ($n$)') 
ax.plot(f_strainzz(sfd_F_ue2), f_J(nlmd_ue2), '--', c=color_n, lw=lw) 

ax.plot(f_strainzz(drex_F_uc1), f_J(drex_nlm_uc1), ':', color=color_b, lw=lw, label=r'D-Rex ($b$)')
ax.plot(f_strainzz(drex_F_ue1), f_J(drex_nlm_ue1), ':', color=color_b, lw=lw)
ax.plot(f_strainzz(drex_F_uc2), f_J(drex_nlm_uc2), ':', color=color_n, lw=lw, label=r'D-Rex ($n$)')
ax.plot(f_strainzz(drex_F_ue2), f_J(drex_nlm_ue2), ':', color=color_n, lw=lw)

ax.plot(np.NaN, np.NaN, '-', color='none', label=' ') # empty entry to order legend rows/cols correctly
ax.plot(np.NaN, np.NaN, '-', color='none', label=' ')

legkwargs = {'handlelength':1.4, 'framealpha':0, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'columnspacing':0.5*1, 'edgecolor':'k', \
            'borderaxespad':0, 'fontsize':FSLEG}
hleg = ax.legend(ncol=2, loc=2, bbox_to_anchor=(0.26,1), **legkwargs)
hleg.get_frame().set_linewidth(0.8)

### Experimental data points

#experiments = (expr_Yabe2020, expr_Boneh2014, expr_Miyazaki2013_uc, expr_Miyazaki2013_ue)
#experiments = (expr_Yabe2020,expr_Miyazaki2013_uc,expr_Miyazaki2013_ue,expr_Boneh2014)
experiments = (expr_Boneh2014,)    
skiplegend = [0,0,1,0]

#mrk = ['s','o','o','^']
mrk = ['d',]
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

h_list = []
leg_list = []
for axname in ['b', 'n']:
    for jj,expr in enumerate(experiments):
    
        x = correlations[jj]['%s20'%(axname)]
        y = correlations[jj]['%s40'%(axname)]
        z = correlations[jj]['%s60'%(axname)]
        J = [f_J(np.array([norm, x[ii],y[ii],z[ii]]), Iend=-1) for ii in np.arange(len(x))]

        X = np.array(expr['Fzz'])
        for ii in np.arange(len(X)):
            J_, x_ = J[ii], X[ii]
            iswarm = expr['T'][ii] > 1250
            lbl = None
            if iswarm: 
                c = color_ln if (axname=='n') else color_lb # warm!
            else:                    
                c = color_dn if (axname=='n') else color_db
                lbl = '%s ($%s$)'%(expr['plotname'], axname)
           
            strainzz = x_ - 1
            h_, = ax.plot(strainzz, J_, ms=ms, markeredgewidth=mew, ls='none', color=c, fillstyle='none', marker=mrk[jj], label=lbl, zorder=1+(len(experiments)-jj))
                
            if ii==0 and not skiplegend[jj] and not iswarm:
                h_list.append(h_)
                leg_list.append(lbl)

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}
leg = plt.legend(h_list, leg_list, fontsize=FSLEG, ncol=1, frameon=False, title_fontsize=FSLEG, \
    loc='lower right', bbox_to_anchor=(1.01,-0.025), bbox_transform=fig.transFigure, **legkwargs); 
leg.get_frame().set_linewidth(0.7);
ax.add_artist(hleg)

### Ticks, limits, etc.

ax.set_xticks(np.arange(-1,15,1))
ax.set_xticks(np.arange(-1,15,0.5), minor=True)
ax.set_yticks(np.arange(1,20,1))
ax.set_yticks(np.arange(1,20,1), minor=True)
ax.set_xlim(xlims)
ax.set_ylim(ylims)

ax.set_ylabel(r'$J$')
ax.set_xlabel(r'$\epsilon_{zz}$')

### Parcel deformation examples

F_uc = sf_F_uc1[np.argmin(np.abs(sf_F_uc1[:,2,2] - 0.5)),:,:]
F_ue = sf_F_ue1[np.argmin(np.abs(sf_F_ue1[:,2,2] - 1.5)),:,:]

W = 0.21
ax1 = plt.axes([0.12, 0.04, W,W], projection='3d')
ax1.patch.set_visible(False)
ax2 = plt.axes([0.33, 0.03, W,W], projection='3d')
ax2.patch.set_visible(False)

kwargs = dict(azim=35, axscale=1, axislabels=True, drawinit=True, fonttex=True, lw=0.75, colorinit=sfplt.c_dred)
sfplt.plotparcel(ax1, F_uc, **kwargs)
sfplt.plotparcel(ax2, F_ue, **kwargs)

### Save fig

f = 'J-index-validation'
plt.savefig('%s.pdf'%(f), dpi=300)
os.system('pdfcrop %s.pdf %s.pdf'%(f,f))

