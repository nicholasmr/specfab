# N. M. Rathmann <rathmann@nbi.ku.dk> 2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmasher as cmr

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=13)
FSSMALL = FS-1
FSBIG = FS+5

#----------------------
# Parameters
#----------------------

nglen = 3
Aglen = 1 # results independent of A because we consider the ratio between power law rheologies

#----------------------
# Initialize model
#----------------------

L = 12
lm, nlm_len = sf.init(L)

m = np.array([0,0,1]) # preferred direction
nlm = sf.nlm_ideal(m, 0, 4) # delta(r-m)

grain_params = sfconst.ice['viscoplastic']['linear']

#----------------------
# Integrate
#----------------------

Na = 100
angles = np.deg2rad(np.linspace(0,90,Na))
Ezz = np.zeros(len(angles))

e1,e2,e3 = np.eye(3)
Eij = sf.Eij_tranisotropic(nlm, e1,e2,e3, *grain_params)
args_rheo = (e1,e2,e3, Eij)

for kk, ang in enumerate(angles):
    Q = R.from_euler('y', ang).as_matrix()
    vw = np.matmul(Q.T,np.matmul(np.tensordot(m, m, axes=0),Q))
    tau = np.eye(3)/3 - vw
    D_G = sf.rheo_fwd_isotropic(  tau, Aglen, nglen)
    D_R = sf.rheo_fwd_orthotropic(tau, Aglen, nglen, *args_rheo)
    Ezz[kk] = np.tensordot(D_R, vw, axes=2) / np.tensordot(D_G, vw, axes=2)

#----------------------
# Plot
#----------------------

di = 3

for ii, angle in enumerate(angles[::di]):
    
    print('theta=%.2f'%(np.rad2deg(angle)))
    
    ### Setup figure

    scale = 0.72
    fig = plt.figure(figsize=(8.1*scale,5*scale))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.85])
    gs.update(hspace=0.0, wspace=0.0, left=0.13, right=1, top=0.98, bottom=0.14)
    ax_E   = fig.add_subplot(gs[0, 0])
    #ax_ODF = fig.add_subplot(gs[0, 1])

    geo, prj = sfplt.getprojection(rotation=-40, inclination=50)
    W = 0.4 
    ax_ODF = plt.axes([0.66, 0.275, W,W], projection=prj) 
    ax_ODF.set_global()

    ylims = [1e-2,2e1]
    X = np.rad2deg(angles)
    lw = 1.5

    ### ODF plot 
    
    cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
    cmap.set_under('w')
    lvlset = [np.linspace(0.15, 1, 7), lambda x,p:'%.1f'%x]
    nlm_rot = sf.rotate_nlm(nlm, angle, 0)
    sfplt.plotODF(nlm_rot, lm, ax_ODF, lvlset=lvlset, cmap=cmap, showcb=False)
    ax_ODF.set_title('ODF', fontsize=FS)

    vi = np.array([ [np.sin(a),0,np.cos(a)] for a in angles[:(ii*di+1)]])
    lat, colat, lon = sfplt.cart2sph(vi, deg=True)
    ax_ODF.plot(lon, lat, '-', lw=2, transform=geo, color=sfplt.c_dgreen)

    kwargs_pnt = dict(ha='center', va='center', transform=geo, color=sfplt.c_dbrown, fontsize=FSBIG)
    lat, colat, lon = sfplt.cart2sph([0,0,1], deg=True)
    ax_ODF.text(lon, lat, r'${\bf \hat{z}}$', **kwargs_pnt)
    sfplt.plotS2point(ax_ODF, [np.sin(angle),0,np.cos(angle)], marker='.', markersize=12, color=sfplt.c_dgreen, transform=geo)
    
    ### Background patches

    ax_E.add_patch(plt.Rectangle((X[0],1), X[-1], ylims[-1], color=sfplt.c_vlblue))
    ax_E.add_patch(plt.Rectangle((X[0],1), X[-1], -1, color=sfplt.c_vlred))

    I = int(len(X)/2-0)
    ax_E.text(X[I], 1.4e0, r'{\bf Softer than isotropy}', c=sfplt.c_dblue, ha='center', va='center', fontsize=FSSMALL-1)
    ax_E.text(X[I], 0.7e0, r'{\bf Harder than isotropy}', c=sfplt.c_dred,  ha='center', va='center', fontsize=FSSMALL-1)

    ax_E.semilogy([np.rad2deg(angle)]*2, ylims, '--', lw=lw,  color=sfplt.c_dgreen)

    ### Plot enahancement factors 

    ax_E.semilogy([0,90], [1,1], ':', lw=lw,  color='k')

    # Shoji and Langway DYE3 deformation tests         
    SL_1890m = np.array([ [0,0.04], [15,0.08], [30,0.29], [90,0.24], [90,0.36],  [45,4.5], [45,6], [45,7] ]).T
    SL_1944m = np.array([ [90,0.24], [90,0.26],   [55,9], [45,6], [45,6.5]]).T
    SL_2006m = np.array([ [0,0.03], [15,0.56], [26,2.8], [45,13], [60,17], [70,3.9], [75, 1.4], [90,0.13] ]).T

    ax_E.semilogy(X, Ezz, '-', color='k', lw=lw, label=r'Orthotropic flow law')

    ms = 6.6
    fc = 'none'
    colors = ['k', 'k', 'k']
    ax_E.semilogy(SL_1890m[0,:], SL_1890m[1,:], '^', fillstyle=fc, markersize=ms, color=colors[2],  label=r'Dye 3, 1890m', clip_on=False)
    ax_E.semilogy(SL_1944m[0,:], SL_1944m[1,:], 's', fillstyle=fc, markersize=ms, color=colors[1],  label=r'Dye 3, 1944m', clip_on=False)
    ax_E.semilogy(SL_2006m[0,:], SL_2006m[1,:], 'o', fillstyle=fc, markersize=ms, color=colors[0],  label=r'Dye 3, 2006m', clip_on=False)

    ax_E.set_xlabel(r'$\theta$')
    ax_E.set_ylabel(r'$\dot{\epsilon}_{zz}/\dot{\epsilon}_{zz}^{\mathrm{Glen}}$')

    xticks = np.arange(0,90+1e-5,5)
    ax_E.set_xticks(xticks[::3])
    ax_E.set_xticks(xticks, minor=True)  
    ax_E.set_xlim([0,90])
    ax_E.set_ylim(ylims)

    legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.75, 'ncol':1, 'columnspacing':0.6, 'handlelength':1.3, 'handletextpad':0.5, 'labelspacing':0.25, 'fontsize':FSSMALL}
    leg = ax_E.legend(loc='lower center', **legkwargs)
    leg.get_frame().set_linewidth(0.8);

    ### Save figure

    fname = 'frames/Emm-misaligned-frame-%i.png'%(ii)
    print('Saving output to %s'%(fname))
    plt.savefig(fname, transparent=1, dpi=250)

