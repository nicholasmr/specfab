#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2023-

import copy, sys, os, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp
import cartopy.crs as ccrs

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=17)

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker

def plot_Ylm(F, lon, lat, ax=None, title='', cmap='RdGy', lvlmax=0.2*4, nchunk=None):
    lvls = np.linspace(-lvlmax,lvlmax,9) # Contour lvls
    hdistr = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='both', nchunk=nchunk, cmap=cmap)
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    ax.set_title(title, pad=10, fontsize=FS)

def discretize(nlm, lm, latres=60):
    theta = np.linspace(0,   np.pi,   latres) # CO-LAT 
    phi   = np.linspace(0, 2*np.pi, 2*latres) # LON
    phi, theta = np.meshgrid(phi, theta) # gridded 
    lon, colat = phi, theta
    lat = np.pi/2-colat
    _,nlm_len = lm.shape
    F = np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], phi,theta) for ii in np.arange(nlm_len) ], axis=0)
    F = np.real(F) # F is real but numerical errors might give small non-zero imag components
    return (F, lon,lat)
    
def newfig(scale=1.0):
    fig = plt.figure(figsize=(1.125*scale,1.37*scale))
    gs = gridspec.GridSpec(1, 1)
    a=0.08
    gs.update(top=0.78, bottom=-0.02, left=a, right=1-a)
    ax = fig.add_subplot(gs[0, 0], projection=prj)
    ax.set_global()
    return ax
    
lw = 1.6
legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.9, 'ncol':1, 'handlelength':1.34, 'labelspacing':0.3}
cbkwargs  = {'orientation':'horizontal', 'fraction':1.3, 'aspect':10}

geo, prj = sfplt.getprojection(rotation=-90-20, inclination=50)
kw_save = dict(transparent=True, dpi=350)

### Expansion

nlm_mag = 1 #1/4
phi = np.deg2rad(0)

lm_list = ((0,0), (2,0), (2,1), (2,2), (4,0), (4,1), (4,2), (4,3), (4,4),)
#lm_list = ((2,2), ) # debug

for (l,m) in lm_list:
    ax = newfig()
    if m == 0:
        (F,lon,lat) = discretize([nlm_mag], np.array([[l,],[m,]]))
        title = '$n_{%i}^{%i} Y_{%i}^{%i}$'%(l,m, l,m)
        plot_Ylm(F, lon, lat, ax=ax, title=title, nchunk=49 if l==0 else None)
    else:
        nlm = nlm_mag*np.exp(m*1j*phi)
        sign = (-1)**(m)
        (F,lon,lat) = discretize([sign*np.conj(nlm),nlm], np.array([[l,l],[m,-m]]))
        title = r'$n_{%i}^{%i}Y_{%i}^{%i} + \mathrm{c.c.}$' % (l,m, l,m)
        plot_Ylm(F, lon, lat, ax=ax, title=title)
    plt.savefig('Ylm/Y%i%i.pdf'%(l,m), **kw_save)        
       
### n(theta,phi) example

ax = newfig()
ax.set_title(r'$n(\theta,\phi)$', pad=10, fontsize=FS)
lm, nlm_len = sf.init(2) 
a2 = np.diag([0.1, 0.3, 0.6])
nlm = sf.a2_to_nlm(a2)
sfplt.plotODF(nlm, lm, ax, lvlset=(np.linspace(0,0.30,9), lambda x,p:'%.1f'%x) , showcb=False, nchunk=None)
plt.savefig('Ylm/nlm.pdf', **kw_save)
