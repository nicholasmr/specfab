# N. M. Rathmann <rathmann@nbi.ku.dk> 2021

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp

sys.path.insert(0, '..')
from header import *
from specfabpy import specfabpy as sf # To use specfabpy compile the specfab Python module by running "make specfabpy"

import warnings
warnings.filterwarnings("ignore")

year2sec = 31556926;

#----------------------
# Parameters
#----------------------

nglen = 3
Aglen = 3.5e-26 # A(T=-25 deg.)
angles = np.deg2rad(np.linspace(0,90,50))
Nt = len(angles)
        
# Linear (n'=1) mixed Taylor--Sachs enhancements: Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
nprime = 1 
Eca   = sf.eca_opt_lin
Ecc   = sf.ecc_opt_lin
alpha = sf.alpha_opt_lin  

#----------------------
# Initialize model
#----------------------

nlm_len  = sf.init(4)         # nlm_len is the number of fabric expansion coefficients (degrees of freedom).
lm       = sf.get_lm(nlm_len) # The (l,m) values corresponding to the coefficients in "nlm".

c_vert = np.array([0,0,1], dtype=np.complex128)
a2_true = np.tensordot(c_vert,  c_vert,  axes=0)
a4_true = np.tensordot(a2_true, a2_true, axes=0)

nlm = sf.a4_to_nlm(a2_true, a4_true)
a2,a4,a6,a8 = sf.ai(nlm)
#print(a2_true)
#print(a2)
print('a2 error: ', np.sum(np.abs(a2_true-a2)))
print('a2 error: ', np.sum(np.abs(a4_true-a4)))

#----------------------
# Integrate
#----------------------

Eang = np.zeros((Nt,3), dtype=np.float64)
e1,e2,e3 = np.array([0,0,1]), np.array([1,0,0]), np.array([0,1,0])

def f_tau(mag, ang):

    c, s = np.cos(ang), np.sin(ang)
    Qy = np.array(((c, 0, s), (0,1,0), (-s,0,c))) # Rotation matrix for rotations around y-axis.

    vw = np.tensordot(c_vert, c_vert, axes=0)    
    vw_rot = np.matmul(np.matmul(Qy, vw), Qy.T)
    
    tau = mag * (np.eye(3)/3 - vw)
    tau_rot = np.matmul(np.matmul(Qy, tau), Qy.T)
    
    return tau_rot, vw_rot

for kk, ang in enumerate(angles):

    Eij = np.transpose(sf.enhfac_eiej(nlm, e1,e2,e3, Ecc, Eca, alpha, nprime))
    tau, vw = f_tau(1, ang)
    eps_G = sf.eps_of_tau__isotropic(         tau, Aglen,nglen)
    eps_R = sf.eps_of_tau__orthotropic(       tau, Aglen,nglen, e1,e2,e3, Eij)
    eps_P = sf.eps_of_tau__orthotropic_pettit(tau, Aglen,nglen, e1,e2,e3, Eij)
    eps_M = sf.eps_of_tau__orthotropic_martin(tau, Aglen,nglen, e1,e2,e3, Eij)
    
    eps_G_vw = np.tensordot(eps_G, vw, axes=2)
    Eang[kk,0] = np.tensordot(eps_R, vw, axes=2)/eps_G_vw
    Eang[kk,1] = np.tensordot(eps_P, vw, axes=2)/eps_G_vw
    Eang[kk,2] = np.tensordot(eps_M, vw, axes=2)/eps_G_vw

#----------------------
# Plot
#----------------------

### Setup figure

scale = 0.69
fig = plt.figure(figsize=(5.2*scale,5.2*scale))
gs = gridspec.GridSpec(1, 1)
ax_E = fig.add_subplot(gs[:, :])

#### Plot enahancement factors 

cglen = '0.4'
ax_E.semilogy([0,90], [1,1], '-', lw=1.2,  color=cglen)
ax_E.text(45, 1.15, "Glen's law", color=cglen, horizontalalignment='center', fontsize=FSSMALL)

# Shoji and Langway DYE3 deformation tests         
SL_1890m = np.array([ [0,0.04], [15,0.08], [30,0.29], [90,0.24], [90,0.36],  [45,4.5], [45,6], [45,7] ]).T
SL_1944m = np.array([ [90,0.24], [90,0.26],   [55,9], [45,6], [45,6.5]]).T
SL_2006m = np.array([ [0,0.03], [15,0.56], [26,2.8], [45,13], [60,17], [70,3.9], [75, 1.4], [90,0.13] ]).T

X = np.rad2deg(angles)
ax_E.semilogy(X, Eang[:,0], '-',  color='k',  label=r'Unapprox.')
ax_E.semilogy(X, Eang[:,2], '--', color='k',  label=r'Martin')
ax_E.semilogy(X, Eang[:,1], ':',  color='k',  label=r'Pettit')
ms = 6.25
fc='none'
#colors = ['#ff7f00', '#1f78b4', '#33a02c']
colors = ['k', 'k', 'k']
ax_E.semilogy(SL_1890m[0,:], SL_1890m[1,:], '^', fillstyle=fc, markersize=ms, color=colors[2],  label=r'DYE 3, 1890m', clip_on=False)
ax_E.semilogy(SL_1944m[0,:], SL_1944m[1,:], 's', fillstyle=fc, markersize=ms, color=colors[1],  label=r'DYE 3, 1944m', clip_on=False)
ax_E.semilogy(SL_2006m[0,:], SL_2006m[1,:], 'o', fillstyle=fc, markersize=ms, color=colors[0],  label=r'DYE 3, 2006m', clip_on=False)
ax_E.set_xlabel(r'$\theta$')
xticks = np.arange(0,90+1e-5,5)
ax_E.set_xticks(xticks[::3])  
ax_E.set_xticks(xticks, minor=True)  
ax_E.set_xlim([0,90])
ax_E.set_ylabel(r'$\dot{\epsilon}_{vv}/\dot{\epsilon}_{vv}^{\mathrm{Glen}}$')
ax_E.set_ylim([1e-2,2e1])
#ax_E.grid()
legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.75, 'ncol':2, 'columnspacing':0.6, 'handlelength':1.3, 'handletextpad':0.5, 'labelspacing':0.25, 'fontsize':FSSMALL-0.5}
leg = ax_E.legend(**legkwargs)
leg.get_frame().set_linewidth(0.8);

### Save figure

fname = 'rheology-compare-Etheta.png'
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=300, bbox_inches = 'tight', pad_inches = 0.05)

