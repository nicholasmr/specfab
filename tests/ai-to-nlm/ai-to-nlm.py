# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2022

"""
Tests conversion from structure-tensor (a2,a4) to spectral expansion coefficients (L=2,L¤4)
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf
            
#----------------------
# Initial, unrotated state
#----------------------

L = 6
lm, nlm_len = sf.init(L)
nlm_0 = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients

# Reference fabric to be rotated
nlm_0[0] = 1/np.sqrt(4*np.pi) # normalized
nlm_0[sf.I20] = 0.5 # n20
nlm_0[sf.I40] = 0.5 # n40
nlm_0[sf.I60] = 0.5 # n60

# Structure tensors of reference fabric
a2_0 = sf.a2(nlm_0)
a4_0 = sf.a4(nlm_0)
a6_0 = sf.a6(nlm_0)

# Rotation
theta = np.deg2rad(45 * 0.8)
phi   = np.deg2rad(45 * 1.2)

#----------------------
# Rotate nlm using Wigner's D matrix to generate nlm with all components non-zero
#----------------------

nlm = sf.rotate_nlm(nlm_0, theta, 0) # vertical rotation
nlm = sf.rotate_nlm(nlm, 0, phi)  # horizontal rotation

a2_wig = sf.a2(nlm)
a4_wig = sf.a4(nlm)
a6_wig = sf.a6(nlm)

print('nlm (Wigner-D rotated): \n', nlm)
#print('a2_wig (tr=%.3f): \n'%(np.trace(a2_wig)), a2_wig)

#----------------------
# Rotate ai using rotation matrices
#----------------------

Rv = np.matrix([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]]) # R_y
Rh = np.matrix([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]]) # R_z

#a2_matrot = np.matmul(Rv,np.matmul(a2_matrot, Rv.T))
#a2_matrot = np.matmul(Rh,np.matmul(a2_matrot, Rh.T))
#print('a2_matrot (tr=%.3f): \n'%(np.trace(a2_matrot)), a2_matrot)

a2 = np.einsum('ki,ij,lj', Rv, a2_0, Rv)
a2 = np.einsum('ki,ij,lj', Rh, a2, Rh)

a4 = np.einsum('mi,nj,ijkl,ok,pl', Rv,Rv, a4_0, Rv,Rv)
a4 = np.einsum('mi,nj,ijkl,ok,pl', Rh,Rh, a4, Rh,Rh)

a6 = np.einsum('oi,pj,qk,ijklmn,rl,sm,tn', Rv,Rv,Rv, a6_0, Rv,Rv,Rv)
a6 = np.einsum('oi,pj,qk,ijklmn,rl,sm,tn', Rh,Rh,Rh, a6, Rh,Rh,Rh)

# Convert to nlm for comparing against wigner-D rotated version 
nlm_a2, nlm_a4, nlm_a6 = np.zeros((nlm_len),dtype=np.complex64), np.zeros((nlm_len),dtype=np.complex64), np.zeros((nlm_len),dtype=np.complex64)
nlm_a2[:sf.L2len] = sf.a2_to_nlm(a2) 
nlm_a4[:sf.L4len] = sf.a4_to_nlm(a4) 
nlm_a6[:sf.L6len] = sf.a6_to_nlm(a6) 

#----------------------
# Print error
#----------------------

np.set_printoptions(linewidth=400)

print('\n--- Conversion error for a2_to_nlm() ---\n')
I = sf.L2len
print(nlm[:I])
print(nlm_a2[:I])
print('\nError (sum abs): ', np.sum(np.abs(nlm[:I]-nlm_a2[:I])))
print('\nDouble check: error (sum abs) of nlm - a2_to_nlm(a2(nlm)): ', np.sum(np.abs(nlm[:I]-sf.a2_to_nlm(sf.a2(nlm))[:I])) )

print('\n--- Conversion error for a4_to_nlm() ---\n')
I = sf.L4len
print(nlm[:I])
print(nlm_a4[:I])
print('\nError (sum abs): ', np.sum(np.abs(nlm[:I]-nlm_a4[:I])))
print('\nDouble check: error (sum abs) of nlm - a4_to_nlm(a4(nlm)): ', np.sum(np.abs(nlm[:I]-sf.a4_to_nlm(sf.a4(nlm))[:I])) )

print('\n--- Conversion error for a6_to_nlm() ---\n')
I = sf.L6len # test full state-vector 
print(nlm[:I])
print(nlm_a6[:I])
print('\nError (sum abs): ', np.sum(np.abs(nlm[:I]-nlm_a6[:I])))
print('\nDouble check: error (sum abs) of nlm - a6_to_nlm(a6(nlm)): ', np.sum(np.abs(nlm[:I]-sf.a6_to_nlm(sf.a6(nlm))[:I])) )

#----------------------
# Plot
#----------------------

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

dpi, scale = 250, 2
fig = plt.figure(figsize=(1.6*scale,2.6*scale))
gs = gridspec.GridSpec(2,1)
a = 0.03
gs.update(left=a, right=1-a/3, top=0.99, wspace=0.015*18, hspace=0.35)
gs.update(wspace=0.015*18, hspace=0.35)

geo = ccrs.Geodetic()
rot0 = 45
rot = 1.0*rot0 # view angle
inclination = 45 # view angle
prj = ccrs.Orthographic(rot, 90-inclination)
ax1 = plt.subplot(gs[0, 0], projection=prj)
ax2 = plt.subplot(gs[1, 0], projection=prj)
ax1.set_global() 
ax2.set_global() 

def plot_axes(ax, geo, cax='tab:red'):
    ax.plot([0],[0],  marker=r'$x$', ms=7, c=cax, transform=geo) # x axis
    ax.plot([90],[0], marker=r'$y$', ms=7, c=cax, transform=geo) # y axis
    ax.plot([0],[90], marker=r'$z$', ms=7, c=cax, transform=geo) # z axis

plot_ODF(nlm, lm, ax=ax1, cmap='Greys', cblabel=r'$\psi/N$ --- Wigner-D rotated nlm; reference case')
plot_axes(ax1, geo)

plot_ODF(nlm_a6, lm, ax=ax2, cmap='Greys', cblabel=r'$\psi/N$ --- a6_to_nlm(R-mat rotated a6)')
plot_axes(ax2, geo)

fout = 'ai-to-nlm.png'
print('\nSaving %s'%(fout))
plt.savefig(fout, dpi=dpi)
    
