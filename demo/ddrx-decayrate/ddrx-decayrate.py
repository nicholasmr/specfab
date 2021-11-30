#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021

import copy, sys, code # code.interact(local=locals())
import numpy as np
sys.path.insert(0, '..')
from header import *
from specfabpy import specfabpy as sf

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.gridspec as gridspec
import scipy.special as sp

latres = 60 # latitude resolution on S^2        
theta = np.linspace(0,   np.pi,   latres) # CO-LAT 
phi   = np.linspace(0, 2*np.pi, 2*latres) # LON
phi, theta = np.meshgrid(phi, theta) # gridded 
lon, colat = phi, theta
lat = np.pi/2-colat

#---------

def plot(clm, ax=None, cmap='RdBu_r', cblbl=r'$\Gamma/\Gamma_0$', titlestr=''):
#    print(np.shape(clm),    clm[1:6],clm[5])
#    clm[1:6] = 0
#    clm[6:] = 0
    ##
    F = np.sum([ clm[ii]*sp.sph_harm(m, l, phi,theta) for ii,(l,m) in enumerate(lm.T) ], axis=0)
    F = np.real(F) # Any imag values are rounding errors.
    llim=0.3
    lvls = np.arange(-llim,llim+1e-5, 0.1)
    hdistr = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), extend='both', cmap=cmap, levels=lvls)
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    cb1 = plt.colorbar(hdistr, ax=ax, fraction=0.055, aspect=10,  orientation='horizontal', pad=0.1, ticks=lvls[1::2])   
    cb1.set_label(cblbl)
    cb1.ax.xaxis.set_ticks(lvls, minor=True)
    ax.set_title(titlestr, fontsize=FS, pad=10)

#---------

inclination = 45 # view angle
rot0 = 1 * -90
rot = 1*-15 + rot0 # view angle

prj = ccrs.Orthographic(rot, 90-inclination)
geo = ccrs.Geodetic()

scale = 2.6
fig = plt.figure(figsize=(3/2*1.45*scale,1.00*scale))
gs = fig.add_gridspec(1,3)
al = 0.04
ar = 0.02
gs.update(left=al, right=1-ar, top=0.98, bottom=0.20, wspace=0.4, hspace=0.4)

ax1 = fig.add_subplot(gs[0,0], projection=prj); ax1.set_global(); 
ax2 = fig.add_subplot(gs[0,1], projection=prj); ax2.set_global(); 
ax3 = fig.add_subplot(gs[0,2], projection=prj); ax3.set_global(); 
axlist = [ax1,ax2,ax3]

#---------

nlm_len = sf.init(4)
lm      = sf.get_lm(nlm_len)
nlm     = np.zeros((nlm_len))
nlm[0]  = 1 # Assume isotropic state for calculating decay rate

### Plot 1
tau = np.diag([0.5,0.5,-1]) # Unconfined UC in z
plot(sf.ddrx_decayrate(nlm, tau), ax1, titlestr=r'%s Unconfined pure shear'%(r'{\Large \bf a}\,\,'))

### Plot 2
tau = np.diag([-1,0,+1]) # Confined UC in z
plot(sf.ddrx_decayrate(nlm, tau), ax2, titlestr=r'%s Confined pure shear'%(r'{\Large \bf b}\,\,'))

### Plot 3
tau = np.matrix([[0,0,1],[0,0,0],[1,0,0]]) # SS xz
#tau = np.matrix([[0,1,0],[1,0,0],[0,0,0]]) # SS xy
plot(sf.ddrx_decayrate(nlm, tau), ax3, titlestr=r'%s Simple shear'%(r'{\Large \bf c}\,\,'))

#---------

for ax in axlist:
    ax.plot([0],[90],'k.', transform=geo)
    ax.plot([rot0],[0],'k.', transform=geo)
    ax.text(rot0-80, 75, r'$\vu{z}$', horizontalalignment='left', transform=geo)
    ax.text(rot0-13, -8, r'$\vu{y}$', horizontalalignment='left', transform=geo)

#---------

fout = 'ddrx-decayrate.png'
print('Saving %s'%(fout))
#plt.tight_layout()
plt.savefig(fout, dpi=300)
#plt.show()

