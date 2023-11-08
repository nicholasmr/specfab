#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

import copy, sys, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()
FSAX = FS+1

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

latres = 60 # latitude resolution on S^2        
colat = np.linspace(0,   np.pi,   latres) 
lon   = np.linspace(0, 2*np.pi, 2*latres) 
lon, colat = np.meshgrid(lon, colat) 
lat = sfdsc.colat2lat(colat)

geo, prj = sfplt.getprojection(rotation=55-90, inclination=50)

def plot(clm, ax=None, cmap='RdBu_r', cblbl=r'$\Gamma/\Gamma_0$', titlestr=''):

    F = np.sum([ clm[ii]*sp.sph_harm(m, l, lon,colat) for ii,(l,m) in enumerate(lm.T) ], axis=0)
    F = np.real(F) # nonzero imag values are rounding errors
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

scale = 2.6
fig = plt.figure(figsize=(3/2*1.45*scale,1.00*scale))
gs = fig.add_gridspec(1,3)
al = 0.04
ar = 0.02
gs.update(left=al, right=1-ar, top=0.98, bottom=0.20, wspace=0.4, hspace=0.4)

ax1 = fig.add_subplot(gs[0,0], projection=prj); ax1.set_global(); 
ax2 = fig.add_subplot(gs[0,1], projection=prj); ax2.set_global(); 
ax3 = fig.add_subplot(gs[0,2], projection=prj); ax3.set_global(); 

lm, nlm_len = sf.init(4) 
nlm = np.zeros((nlm_len))
nlm[0] = 1/np.sqrt(4*np.pi) # Assume isotropic state for calculating decay rate

### Plot 1
tau = np.diag([0.5,0.5,-1]) # Unconfined UC in z
plot(sf.DDRX_decayrate(nlm, tau), ax1, titlestr=r'%s Unconfined pure shear'%(r'{\Large \textit{(a)}}\,\,'))
sfplt.plotcoordaxes(ax1, geo, axislabels='vuxi', fontsize=FSAX, color='k')

### Plot 2
tau = np.diag([0,+1,-1]) # Confined UC in z
plot(sf.DDRX_decayrate(nlm, tau), ax2, titlestr=r'%s Confined pure shear'%(r'{\Large \textit{(b)}}\,\,'))
sfplt.plotcoordaxes(ax2, geo, axislabels='vuxi', fontsize=FSAX, color='k')

### Plot 3
tau = np.einsum('i,j',[0,1,0],[0,0,1]) # y-z shear
tau = (tau + tau.T)/2
plot(sf.DDRX_decayrate(nlm, tau), ax3, titlestr=r'%s Simple shear'%(r'{\Large \textit{(c)}}\,\,'))
sfplt.plotcoordaxes(ax3, geo, axislabels='vuxi', fontsize=FSAX, color='k')

fout = 'ddrx-decayrate.png'
print('Saving %s'%(fout))
plt.savefig(fout, transparent=True, dpi=200)

