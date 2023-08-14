# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
import cmasher as cmr

sys.path.append('../../')
import demolib as dl

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt

FS = sfplt.setfont_tex(fontsize=12)
FSANNO = FS - 1
FSLEG  = FSANNO + 0.5
FSAX   = FSANNO + 2

norm = 1/np.sqrt(4*np.pi) # normalization constant

#--------------------
# Flags
#--------------------

DEBUG = 0 # Low-resolution plotting, etc.

EXPRTYPES = ['uc','ue']
EXPRTYPE = sys.argv[1]
if EXPRTYPE not in EXPRTYPES: raise ValueError('*** Bad experiment type "%s", should be one of: '%(EXPRTYPE),EXPRTYPES) 

#--------------------
# Config
#--------------------

if DEBUG: 
    RESX = RESY = 4*50 
    L = 12
else:
    RESX = RESY = 5*100
    L = 20   

if EXPRTYPE == 'uc': strain_target = -0.90
if EXPRTYPE == 'ue': strain_target = +2

Nt = 100
dn = 4 # save every dn frames (time axis stride)

lm, nlm_len = sf.init(L) 

xlims = [-1.42,2.65]
ylims = [-1.7,3.5]

#--------------------
# Determine enhancements
#--------------------

grain_params = sfconst.ice['viscoplastic']['linear']
(Ezz, Exz, x, y, isvalid) = dl.Eij_statemap_ice(grain_params, xlims, ylims, resx=RESX, resy=RESY)

imdat = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
for ii in range(RESY):
    for jj in range(RESX):
        imdat[ii,jj,:-1] = matplotlib.colors.to_rgb('w') if isvalid[ii,jj] else matplotlib.colors.to_rgb(sfplt.c_lgray)

#--------------------
# Determine correlation line
#--------------------

nl0_unidir, nl0_planar, nl0_circle = sfdsc.nlm_ideal_cases(norm=norm)
x_corr = np.linspace(nl0_planar[0],nl0_unidir[0],100)
y_corr = sf.nhat40_empcorr_ice(x_corr) 

#--------------------
# Load modelled state trajectories
#--------------------

mod = dict(type='pureshear', axis=2, T=1 if strain_target<0 else -1, r=0)
nlm, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=+1, nu=1) # latrot only
nlm /= norm # nlm is now normalized (nhat coefs)
strainzz = np.array([sf.F_to_strain(F[tt,:,:])[2,2] for tt in range(Nt+1)])

print('strain_zz range modelled: ', strainzz[[0,-1]])

#--------------------
# Plot
#--------------------

for tt in np.arange(0,len(strainzz),dn):

    print('%s frame %i, strain_zz = %.3f'%(EXPRTYPE, tt,strainzz[tt]))

    ### Setup figure

    scale = 5.1
    fig = plt.figure(figsize=(1.0*scale,0.7*scale))
    plt.subplots_adjust(left=0.12, right=0.99, top=0.98, bottom=0.16)
    ax = plt.gca()
    cax = fig.add_axes([0.89, 0.19, 0.02, 0.4])
    plt.sca(ax)
        
    ### Plot valid subspace

    im = ax.imshow(imdat, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
    plt.text(1.6, -1.2, '{\\bf Unphysical}\n\\bf{eigenvalues}', color=sfplt.c_dgray, ha='center', fontsize=FSANNO)

    ### Plot Eij

    CS = plt.contourf(x,y, Exz, levels=np.arange(0.5,3.1,0.5), extend='both', cmap=cmr.get_sub_cmap('RdBu', 0.35, 1) , zorder=5)
    cbar = plt.colorbar(CS, cax=cax)
    cbar.ax.set_title('$E_{xz}$', fontsize=FSAX)

    lvls = np.arange(0.2,1.6,0.4)
    CS = plt.contour(x,y,Ezz, levels=lvls, linewidths=0.75, linestyles='-', colors=sfplt.c_vdgray, zorder=10)
    manual_locations = [(0.51, y_) for y_ in np.linspace(-0.5, 1.5, len(lvls))]
    ax.clabel(CS, CS.levels, inline=True, fmt=r'$E_{zz}=%.1f$', fontsize=FSAX, manual=manual_locations)

    ### Plot model state trajectory

    x_, y_ = np.real(nlm[:tt,sf.I20]), np.real(nlm[:tt,sf.I40])
    h_nlm, = ax.plot(x_, y_, ls='-', lw=2.4, c='k', label='Modelled state trajectory', zorder=10)

    ### Plot correlation

    ax.plot(x_corr, y_corr, ':', color=sfplt.c_vdgray, lw=2, label='Empirical correlation', zorder=9)

    ### Misc
    
    dl.plot_nlm_cases(ax, FSANNO, isolbl_above=False)
    
    Eij_grain, alpha, n_grain = grain_params
    plt.text(-1.3, 2.3, r'$(E_{cc}^\prime, E_{ca}^\prime, \alpha) = (%i, 10^{%i}, %.4f)$'%(Eij_grain[0], np.log10(Eij_grain[1]), alpha), color='k', va='center', ha='left', fontsize=FSANNO)

    legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'columnspacing': 0.5, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}
    leg = plt.legend(loc=2, fontsize=FSLEG, frameon=False, ncol=1, **legkwargs); 
    plt.xlabel(r'$\hat{n}_2^0$', fontsize=FSAX)
    plt.ylabel(r'$\hat{n}_4^0$', fontsize=FSAX)
    plt.xlim(xlims)
    plt.ylim(ylims)

    ### Plot ODF inset
    
    if 1:

        if 0: # small ODF inset 
            W = 0.24       
            axpos = [0.1,0.18, W,W]
        else: # big ODF inset
            W = 0.29 # ax width
            axpos = [0.07,0.05, W,W]

        geo, prj = sfplt.getprojection(rotation=50+180, inclination=50)   
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        if EXPRTYPE == 'uc': lvls = np.linspace(0.1,0.9,6)
        if EXPRTYPE == 'ue': lvls = np.linspace(0.1,0.65,6)
        sfplt.plotODF(nlm[tt,:], lm, axin, lvlset=[lvls, lambda x,p:'%.1f'%x], showcb=False);
        axin.set_title('ODF', fontsize=FSLEG)

    ### Save figure

    fname = 'frames/frame-%s-%i.png'%(EXPRTYPE, tt)
    plt.savefig(fname, transparent=1,  dpi=250)
    plt.close()

