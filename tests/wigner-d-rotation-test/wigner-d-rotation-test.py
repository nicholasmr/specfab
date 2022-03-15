# Nicholas M. Rathmann, 2022

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf

### Rotate ODF for L=<4

L = 4
lm, nlm_len = sf.init(L)
nlm = np.zeros((3, nlm_len), dtype=np.complex64) # The expansion coefficients

nlm[0,0] = 1/np.sqrt(4*np.pi) # normalized ODF
nlm[0,3] = 0.25 # vertical single max (nonzero n_2^0 component)

# Wigner D rotation of nlm, implemented in specfab only for components l<=4
lat = np.deg2rad(-45)
lon = np.deg2rad(45)
nlm[1,:] = sf.rotate_nlm4(nlm[0,:], lat, 0) # first rotate around y axis in x--z plane
nlm[2,:] = sf.rotate_nlm4(nlm[1,:], 0, lon) # next rotate around z axis in x--y plane 
            
# Print nlm
np.set_printoptions(precision=4)
I = 6 # l=4 components are identically zero if only n20 is set in the initial (to be rotated) fabric 
print('nlm[:,:I] components up to entries I=%i:'%(I))
print('------------')
print('nlm[0,:I] = ', nlm[0,:I])
print('nlm[1,:I] = ', nlm[1,:I])
print('nlm[2,:I] = ', nlm[2,:I])
print('------------')
print('Im(nlm[1,:I]) = ', np.imag(nlm[1,:I]))
print('Im(nlm[2,:I]) = ', np.imag(nlm[2,:I]))

### Setup figure for plotting

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

def plot_axes(ax, geo):
    cax = 'tab:red'
    ax.plot([0],[0],  marker=r'$x$', ms=7, c=cax, transform=geo) # x axis
    ax.plot([90],[0], marker=r'$y$', ms=7, c=cax, transform=geo) # y axis
    ax.plot([0],[90], marker=r'$z$', ms=7, c=cax, transform=geo) # z axis

dpi, scale = 250, 2.5
fig = plt.figure(figsize=(3*scale,1.4*scale))
gs = gridspec.GridSpec(1,3)
a = 0.03
gs.update(left=a, right=1-a, bottom=0.20, wspace=0.1)

geo = ccrs.Geodetic()
rot = 45 # view angle
inclination = 45 # view angle
prj = ccrs.Orthographic(rot, 90-inclination)
ax1 = plt.subplot(gs[0, 0], projection=prj)
ax2 = plt.subplot(gs[0, 1], projection=prj)
ax3 = plt.subplot(gs[0, 2], projection=prj)

ax1.set_global(); ax2.set_global(); ax3.set_global() 

### Plot

plot_ODF(nlm[0,:], lm, ax=ax1, cmap='Greys', cblabel=r'$\psi/N$')
plot_axes(ax1, geo)

plot_ODF(nlm[1,:], lm, ax=ax2, cmap='Greys', cblabel=r'$\psi/N$')
plot_axes(ax2, geo)

plot_ODF(nlm[2,:], lm, ax=ax3, cmap='Greys', cblabel=r'$\psi/N$')
plot_axes(ax3, geo)

# Save figure
fout = 'wigner-d-rotation-test.png'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=200)
plt.close()

