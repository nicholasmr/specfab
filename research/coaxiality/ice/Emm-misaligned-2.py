# N. M. Rathmann <rathmann@nbi.ku.dk> 2023-2024

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import cmasher as cmr

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc

FS = sfplt.setfont_tex(fontsize=11)
FSSMALL = FS-0

#----------------------
# Parameters
#----------------------

nglen = 3
Aglen = 3.5e-26 # A(T=-25 deg.)

(Eij_grain_0, alpha_0, n_grain_0) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
(Eij_grain_1, alpha_1, n_grain_1) = sfconst.ice['viscoplastic']['linearL10'] # Optimal n'=1 (lin) grain parameters (Rathmann et al., 2024)

print((Eij_grain_1, alpha_1, n_grain_1) )

#----------------------
# Initialize model
#----------------------

L = 10
lm, nlm_len = sf.init(L)

### Unidir state (F_zz = 0)
m = np.array([0,0,1]) # preferred direction
nlm_unidir = sf.nlm_ideal(m, 0, 4) # delta(r-m)

### Uniaxial compression state for L=10 (F_zz = 0.95)
Nt = 200
kwargs_LROT = dict(iota=1, Gamma0=None, nu=1) #      
DK = dict(type='ps', q=+0, tau=+1, axis=2)
strain_target = -0.95 # simulate parcel deformation until this target strain
nlm_uc, F_uc, time_uc, ugrad_uc = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, **kwargs_LROT)
nlm_uc = nlm_uc[-1,:] # pick steady state

#----------------------
# Integrate
#----------------------

angles = np.deg2rad(np.linspace(0,90,150))
Eang = np.zeros((len(angles), 4))

e1,e2,e3 = np.eye(3)
Eij_1 = sf.Eij_tranisotropic(nlm_unidir, e1,e2,e3, Eij_grain_0, alpha_0, n_grain_0)
Eij_2 = sf.Eij_tranisotropic(nlm_uc,     e1,e2,e3, Eij_grain_1, alpha_1, n_grain_1)
Eij_3 = sf.Eij_tranisotropic(nlm_uc,     e1,e2,e3, Eij_grain_1, 0,       n_grain_1)
Eij_4 = sf.Eij_tranisotropic(nlm_uc,     e1,e2,e3, Eij_grain_1, 1,       n_grain_1)

args_rheo_1 = (e1,e2,e3, Eij_1)
args_rheo_2 = (e1,e2,e3, Eij_2)
args_rheo_3 = (e1,e2,e3, Eij_3)
args_rheo_4 = (e1,e2,e3, Eij_4)

def flowlaw_enhancements(tau, vw, nglen):

    D_G = sf.rheo_fwd_isotropic(  tau, Aglen, nglen)
    D_1 = sf.rheo_fwd_orthotropic(tau, Aglen, nglen, *args_rheo_1)
    D_2 = sf.rheo_fwd_orthotropic(tau, Aglen, nglen, *args_rheo_2)
    D_3 = sf.rheo_fwd_orthotropic(tau, Aglen, nglen, *args_rheo_3)
    D_4 = sf.rheo_fwd_orthotropic(tau, Aglen, nglen, *args_rheo_4)

    D_G_vw = np.tensordot(D_G, vw, axes=2)
    Eang = [ np.tensordot(D_, vw, axes=2)/D_G_vw for D_ in (D_1,D_2,D_3,D_4) ]
#    print(Eang)
    return np.array(Eang)

for kk, ang in enumerate(angles):

    Q = R.from_euler('y', ang).as_matrix()
    vw = np.matmul(Q.T,np.matmul(np.tensordot(m, m, axes=0),Q))
    tau = np.eye(3)/3 - vw
    Eang[kk,:] = flowlaw_enhancements(tau, vw, nglen)

#----------------------
# Plot
#----------------------

### Setup figure

scale = 0.69
fig = plt.figure(figsize=(5.2*scale,4.5*scale))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[:, :])

#### Plot enahancement factors 

ylims = [1e-2,2e1]
X = np.rad2deg(angles)
ax.add_patch(plt.Rectangle((X[0],1), X[-1], ylims[-1], color=sfplt.c_vlblue))
ax.add_patch(plt.Rectangle((X[0],1), X[-1], -1, color=sfplt.c_vlred))
ax.text(11, ylims[1]-7, r'$\uparrow$ {\bf Softer}', c=sfplt.c_dblue, ha='center', va='center', fontsize=FSSMALL)
ax.text(11, ylims[0]+0.005, r'$\downarrow$ {\bf Harder}', c=sfplt.c_dred,  ha='center', va='center', fontsize=FSSMALL)

sfplt.panellabel(ax, 1, r'$\epsilon_{zz}=-0.95$', frameon=True, alpha=1.0, fontsize=FS, pad=0.3)

cglen = '0.3'
ax.semilogy([0,90], [1,1], '-', lw=1.1,  color=cglen)
ax.text(45, 1.15-0.53, "Isotropic fabric", color=cglen, horizontalalignment='center', fontsize=FS)

# Shoji and Langway DYE3 deformation tests         
SL_1890m = np.array([ [0,0.04], [15,0.08], [30,0.29], [90,0.24], [90,0.36],  [45,4.5], [45,6], [45,7] ]).T
SL_1944m = np.array([ [90,0.24], [90,0.26],   [55,9], [45,6], [45,6.5]]).T
SL_2006m = np.array([ [0,0.03], [15,0.56], [26,2.8], [45,13], [60,17], [70,3.9], [75, 1.4], [90,0.13] ]).T

X = np.rad2deg(angles)
lw=1.2
#ax.semilogy(X, Eang[:,0], '-',  lw=lw, color='k',          label=r'$\epsilon_{zz}\rightarrow -1$, $\alpha=%.4f$'%(alpha_0))
#ax.semilogy(X, Eang[:,3], ':',  lw=lw, color=sfplt.c_dpurple, label=r'$\epsilon_{zz}=-0.95$, $\alpha=1$')
#ax.semilogy(X, Eang[:,2], '--', lw=lw, color=sfplt.c_dpurple, label=r'$\epsilon_{zz}=-0.95$, $\alpha=0$')
#ax.semilogy(X, Eang[:,1], '-',  lw=lw, color=sfplt.c_dgreen, label=r'$\epsilon_{zz}=-0.95$, $\alpha=%.3f$'%(alpha_1))
ax.semilogy(X, Eang[:,3], ':',  lw=lw*1.1, color=sfplt.c_dpurple, label=r'$\alpha=1$ (Taylor)')
ax.semilogy(X, Eang[:,2], '--', lw=lw, color=sfplt.c_dpurple, label=r'$\alpha=0$ (Sachs)')
ax.semilogy(X, Eang[:,1], '-',  lw=lw, color='k', label=r'$\alpha=%.3f$ (best fit)'%(alpha_1))

ax.fill_between(X, Eang[:,2], y2=Eang[:,3], alpha=0.7, color="#bcbddc")
ax.text(45, 1.95, "Calibration space", color=sfplt.c_dpurple, ma='center', horizontalalignment='center', fontsize=FS)

ms = 6.6
fc = 'none'
colors = ['k', 'k', 'k']
ax.semilogy(SL_1890m[0,:], SL_1890m[1,:], '^', fillstyle=fc, markersize=ms, color=colors[2],  label=r'Dye 3, 1890m', clip_on=False)
ax.semilogy(SL_1944m[0,:], SL_1944m[1,:], 's', fillstyle=fc, markersize=ms, color=colors[1],  label=r'Dye 3, 1944m', clip_on=False)
ax.semilogy(SL_2006m[0,:], SL_2006m[1,:], 'o', fillstyle=fc, markersize=ms, color=colors[0],  label=r'Dye 3, 2006m', clip_on=False)

ax.set_xlabel(r'stress--fabric misalignment $\vartheta$ ($\SI{}{\degree}$)')
#ax.set_ylabel(r'$E_{zz}$')
ax.set_ylabel(r'$\dot{\epsilon}_{zz}/\dot{\epsilon}_{zz}^{(\mathrm{iso})}$')

xticks = np.arange(0,90+1e-5,5)
ax.set_xticks(xticks[::3])  
ax.set_xticks(xticks, minor=True)  
ax.set_xlim([0,90])
ax.set_ylim([1e-2,2e1])

legkwargs = {'frameon':False, 'fancybox':False, 'edgecolor':'k', 'framealpha':1, 'ncol':2, 'columnspacing':0.6, 'handlelength':1.6, 'handletextpad':0.5, 'labelspacing':0.25, 'fontsize':FS}
leg = ax.legend(loc=3, bbox_to_anchor=(-0.05,-0.5), **legkwargs)
leg.get_frame().set_linewidth(0.8);

### Save ODF examples

angles = [30, 75]

for ang in angles:

    geo, prj = sfplt.getprojection(rotation=-30, inclination=45)

    W = H = 0.19

    axpos0 = ang, 3e-2
    trans = ax.transData.transform(axpos0)
    trans = fig.transFigure.inverted().transform(trans)
    axpos = [trans[0]-W/2, trans[1]-H/2, W,H]
    axin = plt.axes(axpos, projection=prj)
    axin.set_global()

    ### Plot ODF    
    lvlset = [np.linspace(0.1, 1.1, 8), lambda x,p:'%.1f'%x]
    cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
    cmap.set_under('w')
    print(np.deg2rad(ang))
    nlm_ = sf.rotate_nlm(nlm_uc, np.deg2rad(ang), 0) 
    sfplt.plotODF(nlm_, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False, nchunk=None)
    
    (e1,e2,e3), _ = sf.eig(nlm_)
    kwargs = dict(fontsize=FSSMALL+1.2, transform=geo)
    sfplt.plotS2text(axin, e1, r'$\vb{m}_1$', c=sfplt.c_green, **kwargs)
    sfplt.plotS2text(axin, [0,0,1], r'$\vu{z}$', c='k', **kwargs)
        
### Save figure

fname = 'Emm-misaligned-2.pdf'
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=175, bbox_inches = 'tight', pad_inches = 0.05)

