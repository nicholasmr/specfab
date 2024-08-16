# Nicholas M. Rathmann, 2022

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

### Rotate ODF for L<=12

L = 12
lm, nlm_len = sf.init(L)
nlm = np.zeros((4, nlm_len), dtype=np.complex64) # The expansion coefficients

a2 = np.diag([.2, .0, .8]) # any second-order structure tensor (not necessarily diagonal)
nlm[0,:sf.L2len] = sf.a2_to_nlm(a2) # l=2 expansion coefficients for corresponding ODF (a2 is normalized)
if 0:
    nlm[0,sf.L4len+2] = 0.2
    nlm[0,sf.L4len+4] = 0.2

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

dpi, scale = 125, 2.0
fig = plt.figure(figsize=(4*scale,1.4*scale))
gs = gridspec.GridSpec(1,4)
a = 0.03
gs.update(left=a, right=1-a, bottom=0.175, wspace=0.1)

geo, prj = sfplt.getprojection(rotation=45, inclination=45)
ax1 = plt.subplot(gs[0, 0], projection=prj)
ax2 = plt.subplot(gs[0, 1], projection=prj)
ax3 = plt.subplot(gs[0, 2], projection=prj)
ax4 = plt.subplot(gs[0, 3], projection=prj)

ax1.set_global(); ax2.set_global(); ax3.set_global(); ax4.set_global()  

### Plot

sfplt.plotODF(nlm[0,:], lm, ax1)
sfplt.plotcoordaxes(ax1, geo, axislabels='vuxi')
ax1.set_title(r'\texttt{nlm}', fontsize=FS)

sfplt.plotODF(nlm[1,:], lm, ax2)
sfplt.plotcoordaxes(ax2, geo, axislabels='vuxi')
ax2.set_title(r'\texttt{nlm\_rot1}', fontsize=FS)

sfplt.plotODF(nlm[2,:], lm, ax3)
sfplt.plotcoordaxes(ax3, geo, axislabels='vuxi')
ax3.set_title(r'\texttt{nlm\_rot2}', fontsize=FS)

sfplt.plotODF(nlm[3,:], lm, ax4)
sfplt.plotcoordaxes(ax4, geo, axislabels='vuxi')
ax4.set_title(r'\texttt{nlm\_rot3}', fontsize=FS)

# Save figure
fout = 'wigner-d-rotation-test.png'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=dpi)
plt.close()

