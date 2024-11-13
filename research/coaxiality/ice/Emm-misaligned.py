# N. M. Rathmann <rathmann@nbi.ku.dk> 2021-2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSSMALL = FS-1

#----------------------
# Parameters
#----------------------

nglen = 3
nglen2 = 4 # for comparison with the recent popular alternative flow exponent
Aglen = 3.5e-26 # A(T=-25 deg.)

(Eij_grain, alpha, n_grain) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)

#----------------------
# Initialize model
#----------------------

lm, nlm_len = sf.init(4)
m = np.array([0,0,1]) # preferred direction
nlm = sf.nlm_ideal(m, 0, 4) # delta(r-m)

#----------------------
# Integrate
#----------------------

angles = np.deg2rad(np.linspace(0,90,50))
Eang = np.zeros((len(angles), 6))

e1,e2,e3 = np.eye(3)
Eij = sf.Eij_tranisotropic(nlm, e1,e2,e3, Eij_grain, alpha, n_grain)
args_rheo = (e1,e2,e3, Eij)

def flowlaw_enhancements(tau, vw, nglen):

    D_G = sf.rheo_fwd_isotropic(         tau, Aglen, nglen)
    D_R = sf.rheo_fwd_orthotropic(       tau, Aglen, nglen, *args_rheo)
    D_P = sf.rheo_fwd_orthotropic_Pettit(tau, Aglen, nglen, *args_rheo)
    D_M = sf.rheo_fwd_orthotropic_Martin(tau, Aglen, nglen, *args_rheo)

    D_G_vw = np.tensordot(D_G, vw, axes=2)
    Eang = [ np.tensordot(D_, vw, axes=2)/D_G_vw for D_ in (D_R,D_P,D_M) ]

    return np.array(Eang)

for kk, ang in enumerate(angles):

    Q = R.from_euler('y', ang).as_matrix()
    vw = np.matmul(Q.T,np.matmul(np.tensordot(m, m, axes=0),Q))
    tau = np.eye(3)/3 - vw
    Eang[kk,:3] = flowlaw_enhancements(tau, vw, nglen)
    Eang[kk,3:] = flowlaw_enhancements(tau, vw, nglen2)

#----------------------
# Plot
#----------------------

### Setup figure

scale = 0.69
fig = plt.figure(figsize=(5.2*scale,5.2*scale))
gs = gridspec.GridSpec(1, 1)
ax_E = fig.add_subplot(gs[:, :])

#### Plot enahancement factors 

cglen = '0.3'
ax_E.semilogy([0,90], [1,1], '-', lw=1.2,  color=cglen)
ax_E.text(45, 1.15, "Glen's law", color=cglen, horizontalalignment='center', fontsize=FSSMALL)

# Shoji and Langway DYE3 deformation tests         
SL_1890m = np.array([ [0,0.04], [15,0.08], [30,0.29], [90,0.24], [90,0.36],  [45,4.5], [45,6], [45,7] ]).T
SL_1944m = np.array([ [90,0.24], [90,0.26],   [55,9], [45,6], [45,6.5]]).T
SL_2006m = np.array([ [0,0.03], [15,0.56], [26,2.8], [45,13], [60,17], [70,3.9], [75, 1.4], [90,0.13] ]).T

X = np.rad2deg(angles)
ax_E.semilogy(X, Eang[:,0], '-',  color='k',    label=r'Unapprox., $n=3$')
ax_E.semilogy(X, Eang[:,3], '-',  color='0.65', label=r'Unapprox., $n=4$')
ax_E.semilogy(X, Eang[:,2], '--', color='k',    label=r'M09, $n=3$')
ax_E.semilogy(X, Eang[:,5], '--', color='0.65', label=r'M09, $n=4$')
ax_E.semilogy(X, Eang[:,1], ':',  color='k',    label=r'P07, any $n$')

ms = 6.6
fc = 'none'
colors = ['k', 'k', 'k']
ax_E.semilogy(SL_1890m[0,:], SL_1890m[1,:], '^', fillstyle=fc, markersize=ms, color=colors[2],  label=r'Dye 3, 1890m', clip_on=False)
ax_E.semilogy(SL_1944m[0,:], SL_1944m[1,:], 's', fillstyle=fc, markersize=ms, color=colors[1],  label=r'Dye 3, 1944m', clip_on=False)
ax_E.semilogy(SL_2006m[0,:], SL_2006m[1,:], 'o', fillstyle=fc, markersize=ms, color=colors[0],  label=r'Dye 3, 2006m', clip_on=False)

ax_E.set_xlabel(r'$\theta$')
ax_E.set_ylabel(r'$\dot{\epsilon}_{vv}/\dot{\epsilon}_{vv}^{\mathrm{Glen}}$')

xticks = np.arange(0,90+1e-5,5)
ax_E.set_xticks(xticks[::3])  
ax_E.set_xticks(xticks, minor=True)  
ax_E.set_xlim([0,90])
ax_E.set_ylim([1e-2,2e1])

legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.75, 'ncol':1, 'columnspacing':0.6, 'handlelength':1.3, 'handletextpad':0.5, 'labelspacing':0.25, 'fontsize':FSSMALL-1.0}
leg = ax_E.legend(**legkwargs, bbox_to_anchor=(0.61,0.54))
leg.get_frame().set_linewidth(0.8);

### Save figure

fname = 'Emm-misaligned.png'
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=175, bbox_inches = 'tight', pad_inches = 0.05)

