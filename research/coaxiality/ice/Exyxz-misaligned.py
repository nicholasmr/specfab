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

FS = sfplt.setfont_tex(fontsize=10)
FSSMALL = FS-0.5

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
mod = dict(type='ps', r=+0, T=+1, axis=0)
strain_target = -0.95 # simulate parcel deformation until this target strain
nlm_uc, F_uc, time_uc, ugrad_uc = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, **kwargs_LROT)
nlm0 = nlm_uc[-1,:] # reference fabric

#----------------------
# Integrate
#----------------------

x,y,z = np.eye(3)
tau1 = (np.einsum('i,j',x,y) + np.einsum('i,j',y,x))
tau2 = -np.sqrt(3)*(np.eye(3)/3 - np.einsum('i,j',y,y))

N = 100
crossover = np.linspace(0,1,N)

D_glen  = np.zeros((N,3,3))
D_orth  = np.zeros((N,3,3))
D_CAFFE = np.zeros((N,3,3))
D_EIE   = np.zeros((N,3,3))

mi  = np.zeros((N,3,3))
nlm = np.zeros((N,nlm_len), dtype=np.complex64)

mi, _ = sfcom.eigenframe(sf.a2(nlm0), modelplane='xy') # are nothing but x,y,z vectors
m1, m2, m3 = mi[:,0], mi[:,1], mi[:,2]
Eij = sf.Eij_tranisotropic(nlm0, m1,m2,m3, *grain_params)
print(m1,m2,m3, Eij)
    
for kk, co in enumerate(crossover):
    tau = (1-co)*tau1 + co*tau2
    D_orth[kk,:,:] = sf.rheo_fwd_orthotropic(tau, Aglen, nglen, m1,m2,m3, Eij)
    D_glen[kk,:,:] = sf.rheo_fwd_isotropic(tau, Aglen, nglen)
    E_CAFFE = sf.E_CAFFE(nlm0, tau, *sfconst.ice['viscoplastic']['CAFFE'])
    E_EIE   = sf.E_EIE(D_orth[kk,:,:], nglen, nlm0, m1,m2,m3, *grain_params)
    D_CAFFE[kk,:,:] = sf.rheo_fwd_isotropic(tau, E_CAFFE*Aglen, nglen)
    D_EIE[kk,:,:]   = sf.rheo_fwd_isotropic(tau,   E_EIE*Aglen, nglen)

### Normalize 
norm = np.sqrt(np.einsum('kij,kji->k', D_glen, D_glen)/2)
print(norm[[0,50,-1]])
def normalize(D): return np.array([D[kk,:,:]/norm[kk] for kk in range(N)])
D_glen, D_orth, D_CAFFE, D_EIE = normalize(D_glen), normalize(D_orth), normalize(D_CAFFE), normalize(D_EIE) 

#----------------------
# Plot
#----------------------

### Setup figure

scale = 0.69
fig = plt.figure(figsize=(4.5*scale,3.0*scale))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[:, :])

#### Plot enahancement factors 

x,y,z = 0,1,2

X = crossover
lw=1.25

c = sfplt.c_dred
kwargs = dict(c=c, lw=lw,)
ax.plot(X[[0,-1]], [Eij[5]]*2, ls=':',  **kwargs)
ax.plot(X[[0,-1]], [Eij[0]]*2, ls=':',  **kwargs)
#ax.plot(X[[0,-1]], [Eij[1]]*2, ls='--', **kwargs)
kwargs = dict(c=c, ha='center', va='center', \
    bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=.21'), clip_on=True, zorder=2)
ax.text(0.8,   Eij[5], r'$E_{12}$', **kwargs)
ax.text(0.4, Eij[0]+0.125, r'$E_{11}$', **kwargs)

kwargs = dict(lw=lw, color='k', zorder=3)
ax.plot(X, D_orth[:,x,y], '-', label=r'$\dot{\epsilon}_{xy}$ Orthotropic', **kwargs)
ax.plot(X, D_orth[:,y,y], '--', label=r'$\dot{\epsilon}_{yy}$ Orthotropic', **kwargs)

kwargs = dict(lw=lw, color='0.52')
ax.plot(X, D_glen[:,x,y], '-',  label=r"$\dot{\epsilon}_{xy}$ Glen", **kwargs)
ax.plot(X, D_glen[:,y,y], '--', label=r"$\dot{\epsilon}_{yy}$ Glen", **kwargs)

cpurp  = '#6a3d9a'
cgreen = 'tab:green'
ax.plot(X, D_CAFFE[:,x,y], '-',  lw=lw, color=cpurp,  label=r"$\dot{\epsilon}_{xy}$ Glen w/ CAFFE")
ax.plot(X, D_CAFFE[:,y,y], '--', lw=lw, color=cpurp,  label=r"$\dot{\epsilon}_{yy}$ Glen w/ CAFFE")
ax.plot(X,   D_EIE[:,x,y], '-',  lw=lw, color=cgreen, label=r"$\dot{\epsilon}_{xy}$ Glen w/ EIE")
ax.plot(X,   D_EIE[:,y,y], '--', lw=lw, color=cgreen, label=r"$\dot{\epsilon}_{yy}$ Glen w/ EIE")

### Remaining axis setup

kwargs = dict(fontsize=FS-0, va='bottom')
y0 = 9.45
ax.text(0, y0, '$\\leftarrow$ $x$--$y$ shear stress', ha='left', **kwargs)
ax.text(1, y0, '$y$-uniaxial tension $\\rightarrow$', ha='right', **kwargs)

ax.annotate('overpredicted \n $y$-extension rate', xy=(0.2, 2.3), xytext=(0.19, 4.6), \
    arrowprops=dict(arrowstyle='-|>, head_width=0.25', facecolor='black'), fontsize=FS-0.6, ha='center', va='center', ma='left')

ax.set_xlabel(r'shear--tension weight', fontsize=FS)
ax.set_ylabel(r'$\dot{\epsilon}_{ij}/\dot{\epsilon}_{\mathrm{iso}}$', fontsize=FS+1)

xticks = np.arange(0,1+1e-5,0.1)
ax.set_xticks(xticks[::2])  
ax.set_xticks(xticks, minor=True)
yticks = np.arange(-5,15,1)
ax.set_yticks(yticks[::1])  
ax.set_yticks(yticks, minor=True)  
ax.set_xlim([0,X[-1]])
ax.set_ylim([-0,9.3])

reorder=lambda hl,nc:(sum((lis[i::nc]for i in range(nc)),[])for lis in hl)
h_l = ax.get_legend_handles_labels()
kwargs_leg = dict(loc=3, bbox_to_anchor=(-0.17,-0.80), frameon=False, fancybox=False, edgecolor='k', framealpha=1, ncol=2, columnspacing=0.8, handlelength=1.1, handletextpad=0.5, labelspacing=0.18, fontsize=FS)
leg = ax.legend(*reorder(h_l, 2), **kwargs_leg)
#leg = ax.legend(**kwargs_leg)

### ODF examples

geo, prj = sfplt.getprojection(rotation=+35, inclination=55)
W = H = 0.23
axpos0 = 0.76, 4.8
trans = ax.transData.transform(axpos0)
trans = fig.transFigure.inverted().transform(trans)
axpos = [trans[0]-W/2, trans[1]-H/2, W,H]
axin = plt.axes(axpos, projection=prj)
axin.set_global()
axin.set_title(r'fixed fabric', fontsize=FSSMALL, pad=5)

lvlset = [np.linspace(0.1, 1.1, 8), lambda x,p:'%.1f'%x]
cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
cmap.set_under('w')
sfplt.plotODF(nlm0, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False, nchunk=None)
sfplt.plotcoordaxes(axin, geo, axislabels='vuxi', negaxes=False, color='k', fontsize=FS) # sfplt.c_dred

kwargs_default = dict(marker='.', ms=6, markeredgewidth=1.0, transform=geo)
ll = 0
sfplt.plotS2point(axin, +mi[:,ll], markerfacecolor=sfplt.c_dred, markeredgecolor=sfplt.c_dred, **kwargs_default)
sfplt.plotS2point(axin, -mi[:,ll], markerfacecolor=sfplt.c_dred, markeredgecolor=sfplt.c_dred, **kwargs_default)
ll = 1
sfplt.plotS2point(axin, +mi[:,ll], markerfacecolor=sfplt.c_dblue, markeredgecolor=sfplt.c_dblue, **kwargs_default)
sfplt.plotS2point(axin, -mi[:,ll], markerfacecolor=sfplt.c_dblue, markeredgecolor=sfplt.c_dblue, **kwargs_default)
#    
kwargs = dict(marker='o', s=8.5, clip_on=False)
X0, Y0, dy = 0.915, 5.3, 0.75
ax.scatter([X0-0.025,],[Y0,],    c=sfplt.c_dred, **kwargs)
ax.scatter([X0-0.025,],[Y0-dy,], c=sfplt.c_dblue, **kwargs)
ax.text(X0, Y0,    r'$\vb{m}_1$', ha='left', va='center')
ax.text(X0, Y0-dy, r'$\vb{m}_2$', ha='left', va='center')
        
### Save figure

fname = 'Exyxz-misaligned.pdf'
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=175, bbox_inches = 'tight', pad_inches = 0.02)

