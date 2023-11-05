#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2023

import numpy as np
import copy, sys, os, code # code.interact(local=locals())
import scipy.special as sp
import cartopy.crs as ccrs

from specfabpy import specfab as sf
from specfabpy import common as sfcom
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=13)
FSLEG = FS-2

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker

#--------------------

def plot_Ylm(F, lon, lat, ax=None, title='', cmap='RdGy', lvlmax=0.2*4):
    lvls = np.linspace(-lvlmax,lvlmax,9) # Contour lvls
    hdistr = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='both', nchunk=5, cmap=cmap)
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    ax.set_title(title, pad=10, fontsize=FS+2)

def discretize_fabric(nlm, lm, latres=60):
    theta = np.linspace(0,   np.pi,   latres) # CO-LAT 
    phi   = np.linspace(0, 2*np.pi, 2*latres) # LON
    phi, theta = np.meshgrid(phi, theta) # gridded 
    lon, colat = phi, theta
    lat = np.pi/2-colat
    _,nlm_len = lm.shape
    F = np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], phi,theta) for ii in np.arange(nlm_len) ], axis=0)
    F = np.real(F) # F is real but numerical errors might give small non-zero imag components
    return (F, lon,lat)
    
lw = 1.6
legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.9, 'ncol':1, 'handlelength':1.34, 'labelspacing':0.3}
cbkwargs  = {'orientation':'horizontal', 'fraction':1.3, 'aspect':10}

### S2 projection
geo, prj = sfplt.getprojection(rotation=-90-20, inclination=50)
 
#-----------------
# n(theta,phi)
#-----------------

scale = 1.3
fig = plt.figure(figsize=(1.125*scale,1.3*scale))
gs = gridspec.GridSpec(1, 1)
a=0.14
gs.update(top=0.80, bottom=0.00, left=a, right=1-a)
ax = fig.add_subplot(gs[0, 0], projection=prj)
ax.set_global()
ax.set_title(r'$n(\theta,\phi)$', pad=10, fontsize=FS+2)

lm, nlm_len = sf.init(2) 
a2 = np.diag([0.1, 0.3, 0.6])
nlm = sf.a2_to_nlm(a2)
#nlm = 1/np.sqrt(4*np.pi) * np.array([1,0.5,-0.325, 0,0,0])
lvlset = (np.linspace(0,0.25,9), lambda x,p:'%.1f'%x) 
sfplt.plotODF(nlm, lm, ax, lvlset=lvlset, showcb=False)
fout = 'nlm.png'
print('Saving %s'%(fout))
plt.savefig(fout, transparent=True, dpi=350)

#-----------------
# Decomposition
#-----------------

nlm_mag = 1 #1/4
phi = np.deg2rad(0)

lm_list = ((0,0), (2,0), (2,1), (2,2), (4,0), (4,1), (4,2), (4,3), (4,4),)
#lm_list = ((2,2), )

for (l,m) in lm_list:

    scale = 1.3
    fig = plt.figure(figsize=(1.125*scale,1.3*scale))
    gs = gridspec.GridSpec(1, 1)
    a=0.14
    gs.update(top=0.80, bottom=0.00, left=a, right=1-a)
    ax = fig.add_subplot(gs[0, 0], projection=prj)
    ax.set_global()
    
    if m == 0:
        (F,lon,lat) = discretize_fabric([nlm_mag], np.array([[l,],[m,]]))
        title = '$n_{%i}^{%i} Y_{%i}^{%i}$'%(l,m, l,m)
        plot_Ylm(F, lon, lat, ax=ax, title=title)

    else:
        nlm = nlm_mag*np.exp(m*1j*phi)
        sign = (-1)**(m)
        (F,lon,lat) = discretize_fabric([sign*np.conj(nlm),nlm], np.array([[l,l],[m,-m]]))
        title = r'$%s\big(n_{%i}^{%i}\big)^* Y_{%i}^{-%i} + n_{%i}^{%i}Y_{%i}^{%i}$' % ('+' if sign>0 else '-', l,m, l,m, l,m, l,m)
        plot_Ylm(F, lon, lat, ax=ax, title=title)

    fout = 'Y%i%i.png'%(l,m)
    print('Saving %s'%(fout))
    plt.savefig(fout, transparent=True, dpi=350)
    
   
