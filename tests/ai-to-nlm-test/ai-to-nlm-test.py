# Nicholas M. Rathmann, 2022

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf
            
#----------------------
# Initial, unrotated state
#----------------------

lm, nlm_len = sf.init(4)
nlm_0 = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients

# Reference fabric to be rotated
nlm_0[0]  = 1/np.sqrt(4*np.pi) # normalized
nlm_0[3]  = 0.5 # n20
nlm_0[10] = 0.5 # n40

# Structure tensors of reference fabric
a2_0 = sf.a2(nlm_0)
a4_0 = sf.a4(nlm_0)

# Rotation
theta = np.deg2rad(45 * 1)
phi   = np.deg2rad(45 * 1)

#----------------------
# Rotate nlm using Wigner's D matrix
#----------------------

nlm_wig = sf.rotate_nlm4(nlm_0, -theta, 0) # vert
nlm_wig = sf.rotate_nlm4(nlm_wig, 0, phi) # horiz
a2_wig = sf.a2(nlm_wig)
a4_wig = sf.a4(nlm_wig)

print('a2_wig (tr=%.3f): \n'%(np.trace(a2_wig)), a2_wig)

#----------------------
# Rotate ai using rotation matrices
#----------------------

a2_matrot = a2_0
a4_matrot = a4_0

Rv = np.matrix([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]]) # R_y
#Rv = np.matrix([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]]) # R_x not used
Rh = np.matrix([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]]) # R_z

a2_matrot = np.matmul(Rv,np.matmul(a2_matrot, Rv.T))
a2_matrot = np.matmul(Rh,np.matmul(a2_matrot, Rh.T))
print('a2_matrot (tr=%.3f): \n'%(np.trace(a2_matrot)), a2_matrot)

a4_matrot = np.einsum('mi,nj,ok,pl,ijkl', Rv, Rv, Rv, Rv, a4_matrot)
a4_matrot = np.einsum('mi,nj,ok,pl,ijkl', Rh, Rh, Rh, Rh, a4_matrot)

nlm_matrot = np.zeros((nlm_len), dtype=np.complex64) 
#nlm_matrot[:6] = sf.a2_to_nlm(a2_matrot) # uncomment to instead test a2_to_nlm()
nlm_matrot[:] = sf.a4_to_nlm(a4_matrot)

#----------------------
# Print error
#----------------------

np.set_printoptions(linewidth=400)
print('--- Error ---')
print(nlm_wig[:])
print(nlm_matrot[:])
print('nlm err (sum abs): ', np.sum(np.abs(nlm_wig-nlm_matrot)))

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

def plot_axes(ax, geo, cax='#ff7f00'):
    ax.plot([0],[0],  marker=r'$x$', ms=7, c=cax, transform=geo) # x axis
    ax.plot([90],[0], marker=r'$y$', ms=7, c=cax, transform=geo) # y axis
    ax.plot([0],[90], marker=r'$z$', ms=7, c=cax, transform=geo) # z axis

plot_ODF(nlm_wig, lm, ax=ax1, cmap='Greys', cblabel=r'$\psi/N$ (rot. nlm)')
plot_axes(ax1, geo)

plot_ODF(nlm_matrot, lm, ax=ax2, cmap='Greys', cblabel=r'$\psi/N$ (rot. ai)')
plot_axes(ax2, geo)

fout = 'ai-to-nlm-test.png'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=dpi)
    
