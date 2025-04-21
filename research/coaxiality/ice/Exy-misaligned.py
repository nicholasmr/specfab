# Nicholas Rathmann <rathmann@nbi.ku.dk> 2024

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from scipy.spatial.transform import Rotation as R

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import cmasher as cmr

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import common as sfcom
from specfabpy import plotting as sfplt
from specfabpy import integrator as sfint

FS = sfplt.setfont_tex(fontsize=9.5)
FSSMALL = FS

#----------------------
# Parameters
#----------------------

nglen = 3
Aglen = 3.5e-26 # A(T=-25 deg.)
grain_params = sfconst.ice['viscoplastic']['linearL10'] # Optimal n'=1 (lin) grain parameters (Rathmann et al., 2024)

#----------------------
# Initialize model
#----------------------

L = 10
lm, nlm_len = sf.init(L)

# Uniaxial compression fabric as reference fabric state
Nt = 200
kwargs_LROT = dict(iota=1, Gamma0=None, nu=1) #      
DK = dict(type='ps', q=+0, tau=+1, axis=1)
strain_target = -0.95 # simulate parcel deformation until this target strain
nlm_uc, F_uc, time_uc, ugrad_uc = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, **kwargs_LROT)
nlm0 = nlm_uc[-1,:] # reference fabric

#----------------------
# Integrate
#----------------------

x,y,z = np.eye(3)
tau = 1*(np.einsum('i,j',x,y) + np.einsum('i,j',y,x))

N = 100
anglesdeg = np.linspace(0,45,N)
anglesrad = np.deg2rad(anglesdeg)

epsij_glen = sf.rheo_fwd_isotropic(tau, Aglen, nglen)
epsij_orth = np.zeros((N,3,3))
norm = np.sqrt(np.einsum('ij,ji',epsij_glen,epsij_glen)/2) # ||eps_ij^(Glen)||
mi  = np.zeros((N,3,3))
nlm = np.zeros((N,nlm_len), dtype=np.complex64)
E_CAFFE = np.zeros((N))
E_EIE   = np.zeros((N))

for kk, ang in enumerate(anglesrad):
    nlm[kk,:] = sf.rotate_nlm(nlm0, 0, ang)
    nlm_ = nlm[kk,:].copy()
    mi[kk,:,:], _ = sfcom.eigenframe(sf.a2(nlm_), modelplane='xy') # are nothing but x,y,z vectors
    Eij = sf.Eij_tranisotropic(nlm_, *mi[kk], *grain_params)
    #print(ang, mi[kk], Eij)
    epsij_orth[kk,:,:] = sf.rheo_fwd_orthotropic(tau, Aglen, nglen, *mi[kk], Eij)
    E_CAFFE[kk] = sf.E_CAFFE(nlm_, tau, *sfconst.ice['viscoplastic']['CAFFE'])
    E_EIE[kk]   = sf.E_EIE(epsij_orth[kk,:,:], nglen, nlm_, *mi[kk], *grain_params)

#----------------------
# Plot
#----------------------

### Setup figure

scale = 0.69
fig = plt.figure(figsize=(4.8*scale,3.0*scale))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[:, :])

#### Plot enahancement factors 

x,y,z = 0,1,2

X = anglesdeg
lw=1.25

Z = epsij_orth/norm
kwargs = dict(lw=lw, color='k')
ax.plot(X, Z[:,x,y], '-', label=r'$\dot{\epsilon}_{xy}$ Orthotropic', **kwargs)
ax.plot(X, Z[:,y,y], '--', label=r'$\dot{\epsilon}_{yy}$ Orthotropic', **kwargs)

Z = epsij_glen/norm
ax.plot(X, [Z[x,y]]*len(X), '-',  lw=lw, color='0.52', label=r"$\dot{\epsilon}_{xy}$ Glen", zorder=3)
ax.plot(X, [Z[y,y]]*len(X), '--',  lw=lw, color='0.52', label=r"$\dot{\epsilon}_{yy}$ Glen", zorder=3, clip_on=False)

c = sfplt.c_dred
ax.plot(X[[0,-1]], [Eij[5]]*2, c=c, lw=lw, ls=':')
ax.plot(X[[0,-1]], [Eij[0]]*2, c=c, lw=lw, ls=':')
#ax.plot(X[[0,-1]], [Eij[1]]*2, c=c, lw=lw, ls='--')

kwargs = dict(c=c, ha='center', va='center', \
    bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=.21'), clip_on=True, zorder=2)
ax.text(16,   Eij[5], r'$E_{12}$', **kwargs)
ax.text(16, Eij[0]+0.125, r'$E_{11}$', **kwargs)

cpurp  = '#6a3d9a'
cgreen = 'tab:green'
ax.plot(np.NaN, np.NaN, '-', color='none', label='')
ax.plot(np.NaN, np.NaN, '-', color='none', label='')
ax.plot(X, E_CAFFE, '-',  lw=lw, color=cpurp, label=r"$\dot{\epsilon}_{xy}$ Glen w/ CAFFE")
ax.plot([0,1], [-10]*2, '', color='none', label=r" ") # dummy to adjust legend entries
ax.plot(X, E_EIE, '-',  lw=lw, color=cgreen, label=r"$\dot{\epsilon}_{xy}$ Glen w/ EIE")

### Remaining axis setup

ax.text(32, 7.3, 'constant $x$--$y$ shear stress', fontsize=FSSMALL, ma='center', ha='center', va='center')

kwargs = dict(fontsize=FSSMALL, va='bottom')
y0 = 9.42
ax.text(0, y0, '$\\leftarrow$ compatible', ha='left', **kwargs)
ax.text(45, y0, 'incompatible $\\rightarrow$', ha='right', **kwargs)

ax.annotate('noncoaxial\n behaviour', xy=(12, 2.1), xytext=(7, 4.2), \
    arrowprops=dict(arrowstyle='-|>, head_width=0.25', facecolor='black'), fontsize=FS, ha='center', va='center', ma='center')

ax.set_xlabel(r'stress--fabric misalignment angle $\vartheta$ ($\SI{}{\degree}$)', fontsize=FS)
ax.set_ylabel(r'$\dot{\epsilon}_{ij}/\dot{\epsilon}_{\mathrm{iso}}$', fontsize=FS+0.5)

xticks = np.arange(0,90+1e-5,5)
ax.set_xticks(xticks[::2])  
ax.set_xticks(xticks, minor=True)
yticks = np.arange(-5,15,1)
ax.set_yticks(yticks[::1])  
ax.set_yticks(yticks, minor=True)  
ax.set_xlim([0,X[-1]])
ax.set_ylim([-0.01,9.3])

reorder=lambda hl,nc:(sum((lis[i::nc]for i in range(nc)),[])for lis in hl)
h_l = ax.get_legend_handles_labels()
kwargs_leg = dict(loc=3, bbox_to_anchor=(-0.05,-0.78), frameon=False, fancybox=False, edgecolor='k', framealpha=1, ncol=2, columnspacing=0.8, handlelength=1.1, handletextpad=0.5, labelspacing=0.18, fontsize=FS)
leg = ax.legend(*reorder(h_l, 2), **kwargs_leg)
#leg = ax.legend(**kwargs_leg)


### ODF examples

angles_plot = [0,20,40]
I = [np.argmin(np.abs(a-anglesdeg)) for a in angles_plot]

for kk in I:
    
    ang = anglesdeg[kk]
    geo, prj = sfplt.getprojection(rotation=+35-90, inclination=55)
    W = H = 0.23
    axpos0 = ang, 11.9
    trans = ax.transData.transform(axpos0)
    trans = fig.transFigure.inverted().transform(trans)
    axpos = [trans[0]-W/2, trans[1]-H/2, W,H]
    axin = plt.axes(axpos, projection=prj)
    axin.set_global()
    axin.set_title(r'$\vartheta=\SI{%i}{\degree}$'%(ang), fontsize=FS, pad=3)

    lvlset = [np.linspace(0.1, 1.1, 8), lambda x,p:'%.1f'%x]
    cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
    cmap.set_under('w')
    sfplt.plotODF(nlm[kk,:], lm, axin, lvlset=lvlset, cmap=cmap, showcb=False, nchunk=None)
    sfplt.plotcoordaxes(axin, geo, axislabels='vuxi', negaxes=True, color='k', fontsize=FS) # sfplt.c_dred

    kwargs_default = dict(marker='.', ms=6, markeredgewidth=1.0, transform=geo)
    ll = 0
    sfplt.plotS2point(axin, +mi[kk,ll], markerfacecolor=sfplt.c_dred, markeredgecolor=sfplt.c_dred, **kwargs_default)
    sfplt.plotS2point(axin, -mi[kk,ll], markerfacecolor=sfplt.c_dred, markeredgecolor=sfplt.c_dred, **kwargs_default)
    ll = 1
    sfplt.plotS2point(axin, +mi[kk,ll], markerfacecolor=sfplt.c_dblue, markeredgecolor=sfplt.c_dblue, **kwargs_default)
    sfplt.plotS2point(axin, -mi[kk,ll], markerfacecolor=sfplt.c_dblue, markeredgecolor=sfplt.c_dblue, **kwargs_default)
    
kwargs = dict(marker='o', s=8.5, clip_on=False)
X0, Y0, dy = 7.5, 12.6, 0.75
ax.scatter([X0-1.4,],[Y0,],    c=sfplt.c_dred, **kwargs)
ax.scatter([X0-1.4,],[Y0-dy,], c=sfplt.c_dblue, **kwargs)
ax.text(X0, Y0,    r'$\vb{m}_1$', ha='left', va='center')
ax.text(X0, Y0-dy, r'$\vb{m}_2$', ha='left', va='center')
        
### Save figure

fname = 'Exy-misaligned.pdf'
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=175, bbox_inches = 'tight', pad_inches = 0.02)

