# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
from progress.bar import Bar

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmasher as cmr

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

normfac = 1/np.sqrt(4*np.pi) # normalization constant

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

xlims, ylims = [-1.42,2.65], [-1.7,3.5]
sc = np.diff(ylims)/np.diff(xlims)
x = np.linspace(xlims[0],xlims[1],RESX)
y = np.linspace(ylims[0],ylims[1],RESY)

#--------------------
# Determine ideal boundary line
#--------------------

m = [0,0,1]
Il24 = [sf.I20, sf.I40] # l=2,4, m=0 coefs
n20_unidir, n40_unidir = np.real(sf.nlm_ideal(m, 0, L))[Il24]/normfac       # delta distributed
n20_planar, n40_planar = np.real(sf.nlm_ideal(m, np.pi/2, L))[Il24]/normfac # x--y planar distributed
n20_circle, n40_circle = np.real(sf.nlm_ideal(m, np.pi/4, L))[Il24]/normfac # 45 deg circle distributed

#--------------------
# Determine valid subspace
#--------------------

xv, yv = np.meshgrid(x, y, indexing='xy')
isvalid = np.reshape(sf.nlm_isvalid(xv.flatten(), yv.flatten()), (RESY, RESX))
imdat = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
for ii in range(RESY):
    for jj in range(RESX):
        imdat[ii,jj,:-1] = [1]*3 if isvalid[ii,jj] else [0.85]*3

#--------------------
# Determine enhancements
#--------------------

Ezz = np.zeros((RESY, RESX)) 
Exz = np.zeros((RESY, RESX)) 
(Eij_grain, alpha, n_grain) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)

m, t = np.array([0,0,1]), np.array([1,0,0])
mm, mt = np.tensordot(m,m, axes=0), np.tensordot(m,t, axes=0)
tau_mm = 1*(np.identity(3)-3*mm)
tau_mt = mt + mt.T

for ii, y_ in enumerate(y):
    for jj, x_ in enumerate(x): 
        if isvalid[ii,jj]:
            nlm_ = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients
            nlm_[0], nlm_[sf.I20], nlm_[sf.I40] = 1, x_, y_
            Ezz[ii,jj] = sf.Evw_tranisotropic(nlm_, m,m,tau_mm, Eij_grain,alpha,n_grain)
            Exz[ii,jj] = sf.Evw_tranisotropic(nlm_, m,t,tau_mt, Eij_grain,alpha,n_grain)
        else:
            Ezz[ii,jj] = np.nan
            Exz[ii,jj] = np.nan

#--------------------
# Determine correlation line
#--------------------

x_corr = np.linspace(n20_planar,n20_unidir,100)
y_corr = sf.nhat40_empcorr_ice(x_corr) 

#--------------------
# Load modelled state trajectories
#--------------------

mod = dict(type='pureshear', axis=2, T=1 if strain_target<0 else -1, r=0)
nlm, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=+1, nu=1) # latrot only
nlm /= normfac # nlm is now normalized (nhat coefs)
strainzz = np.array([sf.F_to_strain(F[tt,:,:])[2,2] for tt in range(Nt+1)])

print('strain_zz range modelled: ', strainzz[[0,-1]])

#--------------------
# Plot
#--------------------

ms = 7
mse = ms+0.5 # end-member case points
FSANNO = FS - 1
FSLEG = FSANNO + 0.5
FSAXLBL = FSANNO + 2
FSEij  = FSANNO + 2

cdarkgray = '0.125'

c_planar = '#8c510a' 
c_unidir = '#01665e'
c_circle = '#762a83'

cl_unidir = np.array([199,234,229])/255
cl_planar = np.array([246,232,195])/255

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'columnspacing': 0.5, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}


for tt in np.arange(0,len(strainzz),dn):

    print('%s frame %i, strain_zz = %.3f'%(EXPRTYPE, tt,strainzz[tt]))

    ### Setup figure

    scale = 5.1
    fig = plt.figure(figsize=(1.0*scale,0.7*scale))
    plt.subplots_adjust(left=0.12, right=0.99, top=0.98, bottom=0.16)
    ax = plt.gca()

    ### Plot valid subspace

    im = ax.imshow(imdat, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
    x_,y_ = -0.99, 1.9
    x_,y_ = 1.6, -1.2
    plt.text(x_,y_, '{\\bf Unphysical}\n\\bf{eigenvalues}', color='0.3', ha='center', rotation=0, fontsize=FSANNO)

    ### Plot Eij

    cmap = cmr.get_sub_cmap('RdBu', 0.35, 1) 
    CS = plt.contourf(x,y, Exz, levels=np.arange(0.5,3.1,0.5), extend='both', cmap=cmap, zorder=5)
    cax = fig.add_axes([0.89, 0.19, 0.02, 0.4])
    cbar = plt.colorbar(CS, cax=cax)
#    cbar.ax.set_ylabel('$E_{xz}$')
    cbar.ax.set_title('$E_{xz}$', fontsize=FSEij)
    plt.sca(ax)

    lvls = np.arange(0.2,1.6,0.4)
    CS = plt.contour(x,y,Ezz, levels=lvls, linewidths=0.75, linestyles='-', colors=cdarkgray, zorder=10)
    manual_locations = [(0.51, y_) for y_ in np.linspace(-0.5, 1.5, len(lvls))]
    ax.clabel(CS, CS.levels, inline=True, fmt=r'$E_{zz}=%.1f$', fontsize=FSEij, manual=manual_locations)

    ### End-member cases

    dy = 0.04/normfac

    # Isotropic state
    ax.plot(0,0,'o', ms=mse, c=cdarkgray, zorder=20)
    plt.text(0, 0-dy, r'{\bf Isotropic}', color=cdarkgray, ha='center', va='top', fontsize=FSANNO, zorder=20)

    # Unidirectional
    ax.plot(n20_unidir,n40_unidir, marker='o', ms=mse, ls='none', c=c_unidir, zorder=20)
    plt.text(n20_unidir-0.1, n40_unidir+dy, '{\\bf Unidirectional}', color=c_unidir, ha='center', va='bottom', ma='center', fontsize=FSANNO, zorder=20)

    # Planar
    ax.plot(n20_planar, n40_planar, marker='o', ms=mse, ls='none', c=c_planar, zorder=20)
    plt.text(n20_planar, n40_planar+dy, '{\\bf Planar}', color=c_planar, ha='center', va='bottom', ma='center', fontsize=FSANNO, zorder=20)

    # Circle
    ax.plot(n20_circle, n40_circle, marker='o', ms=mse, ls='none', c=c_circle, zorder=20)
    plt.text(n20_circle, n40_circle-dy, '{\\bf Circle}', color=c_circle, ha='center', va='top', ma='center', fontsize=FSANNO, zorder=20)

    ### Plot model state trajectory

    x_, y_ = np.real(nlm[:tt,sf.I20]), np.real(nlm[:tt,sf.I40])
    h_nlm, = ax.plot(x_, y_, ls='-', lw=2.4, c='k', label='Modelled state trajectory', zorder=10)

    ### Plot correlation

    ax.plot(x_corr, y_corr, ':', color=cdarkgray, lw=2, label='Empirical correlation', zorder=9)

    ### Aux
    
    plt.text(-1.3, 2.3, r'$(E_{cc}^\prime, E_{ca}^\prime, \alpha) = (%i, 10^{%i}, %.4f)$'%(Eij_grain[0], np.log10(Eij_grain[1]), alpha), color='k', va='center', ha='left', fontsize=FSANNO)

    leg = plt.legend(loc=2, fontsize=FSLEG, frameon=False, ncol=1, **legkwargs); 
    plt.sca(ax)
    plt.xlabel(r'$\hat{n}_2^0$', fontsize=FSAXLBL)
    plt.ylabel(r'$\hat{n}_4^0$', fontsize=FSAXLBL)
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

