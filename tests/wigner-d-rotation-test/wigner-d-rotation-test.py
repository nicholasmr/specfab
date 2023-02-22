# Nicholas M. Rathmann, 2022

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf

### Rotate ODF for L<=12

L = 12
lm, nlm_len = sf.init(L)
nlm = np.zeros((4, nlm_len), dtype=np.complex64) # The expansion coefficients

if 1:
    a2 = np.diag([0.0,0.0,1.0]) # any second-order structure tensor (not necessarily diagonal)
    nlm[0,0:6] = sf.a2_to_nlm(a2) # l=2 expansion coefficients for corresponding ODF (a2 is normalized)
else:
    nlm[0,0] = 1/np.sqrt(4*np.pi)
#    nlm[0,nlm_len-13] = 0.25
    nlm[0,nlm_len-(2*12+1)-(2*10+1)-0*(2*8+1)-0*(2*6+1)-0*(2*4+1)] = 0.3

# Wigner D rotation of nlm, implemented in specfab only for components l<=12
lat = np.deg2rad(60) # z-->x angle
lon = np.deg2rad(45) # x-->y angle
nlm[1,:] = sf.rotate_nlm(nlm[0,:], lat, 0) # first rotate around y axis in x--z plane
nlm[2,:] = sf.rotate_nlm(nlm[1,:], 0, lon) # next rotate around z axis in x--y plane 

nlm[3,:] = sf.rotate_nlm(nlm[2,:], -lat, -lon) # rotate back
            
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

dpi, scale = 125, 2.2
fig = plt.figure(figsize=(4*scale,1.4*scale))
gs = gridspec.GridSpec(1,4)
a = 0.03
gs.update(left=a, right=1-a, bottom=0.175, wspace=0.1)

geo = ccrs.Geodetic()
rot = 45 # view angle
inclination = 45 # view angle
prj = ccrs.Orthographic(rot, 90-inclination)
ax1 = plt.subplot(gs[0, 0], projection=prj)
ax2 = plt.subplot(gs[0, 1], projection=prj)
ax3 = plt.subplot(gs[0, 2], projection=prj)
ax4 = plt.subplot(gs[0, 3], projection=prj)

ax1.set_global(); ax2.set_global(); ax3.set_global(); ax4.set_global()  

### Plot

plot_ODF(nlm[0,:], lm, ax=ax1, cmap='Greys', cblabel=r'$\psi/N$')
plot_axes(ax1, geo)
ax1.set_title('nlm')

plot_ODF(nlm[1,:], lm, ax=ax2, cmap='Greys', cblabel=r'$\psi/N$')
plot_axes(ax2, geo)
ax2.set_title('nlm_rot1')

plot_ODF(nlm[2,:], lm, ax=ax3, cmap='Greys', cblabel=r'$\psi/N$')
plot_axes(ax3, geo)
ax3.set_title('nlm_rot2')

plot_ODF(nlm[3,:], lm, ax=ax4, cmap='Greys', cblabel=r'$\psi/N$')
plot_axes(ax4, geo)
ax4.set_title('nlm_rot3 (back)')

# Save figure
fout = 'wigner-d-rotation-test.png'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=dpi)
plt.close()

