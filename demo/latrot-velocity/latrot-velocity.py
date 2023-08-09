#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

import copy, sys, code # code.interact(local=locals())
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import cartopy.crs as ccrs

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSSMALL = FS-1

### High-res grid

mul = 10
colat = np.linspace(0, np.pi, 20*mul)   
lon   = np.linspace(0, 2*np.pi, 10*mul)
lon2, colat2 = np.meshgrid(lon, colat)
rv, tv, pv = sfdsc.sphericalbasisvectors(colat2, lon2)

### Low-res grid

mul = 1.0
colat_ = np.linspace(0, np.pi, int(10*mul))   
lon_   = np.linspace(0, 2*np.pi, int(20*mul))
lon2_, colat2_ = np.meshgrid(lon_, colat_)
rv_, tv_, pv_ = sfdsc.sphericalbasisvectors(colat2_, lon2_)

### Functions

def get_velmap(ugrad, rv, tv, pv, speedthres=0):
    D, W = (ugrad + ugrad.T)/2, (ugrad - ugrad.T)/2
    Wp = np.einsum('ijk,ljk,lm,mjk->ijk', rv,rv,D,rv) - np.einsum('ij,jkl->ikl', D, rv)
    cdot = np.einsum('ij,jkl->ikl', W, rv) + Wp
    cdot *= -1 # @TODO there is a sign error somewhere, but fixed here (to be found)
    ut = np.einsum('ijk,ijk->jk', cdot, tv)
    up = np.einsum('ijk,ijk->jk', cdot, pv)
    speed = np.linalg.norm(cdot, axis=0)
    ut[speed<speedthres] = np.nan
    return (ut, up, speed)
    

def plot(ugrad, ax, titlestr='', speedthres=0):

    transform = ccrs.PlateCarree()

    ### Contour, high-res

    ut, up, speed = get_velmap(ugrad, rv, tv, pv, speedthres=speedthres)
    lvls = np.linspace(0,0.8,5)
    F = speed/np.linalg.norm(ugrad)
    x, y = np.rad2deg(lon), np.rad2deg(sfdsc.colat2lat(colat))
    hdistr = ax.contourf(x, y, F, transform=transform, extend='max', cmap='YlOrBr', levels=lvls)

    ### Quiver, low-res

    ut_, up_, speed_ = get_velmap(ugrad, rv_, tv_, pv_, speedthres=speedthres)
    lvls = np.linspace(0,0.8,5)
    x_, y_ = np.rad2deg(lon_), np.rad2deg(sfdsc.colat2lat(colat_)) # np.rad2deg(colat_)
    if 0: # debug/verify quiver 
        x_, y_ = np.array([180]*3), sfdsc.colat2lat(np.array([0, 45, 90]), deg=True)
        up_, ut_ = np.array([1]*3), np.array([0]*3)
    QV1 = ax.quiver(x_,y_, up_,ut_, scale=7.5, width=0.012, color='k', transform=transform)
    plt.quiverkey(QV1, 0.95, 0.00, 1, r'$\dot{\vb{c}}$', labelpos='E', coordinates='axes', labelsep=0.05)

    if 0: # debug cordinate transform
        ax.plot(180, 90- 0, 'Xr', transform=transform)
        ax.plot(180, 90- 45, 'or', transform=transform)
        ax.plot(180, 90- 90, 'sr', transform=transform)
        ax.plot(180*0, 90- 0, 'Xb', transform=transform)
        ax.plot(180*0, 90- 45, 'ob', transform=transform)
        ax.plot(180*0, 90- 90, 'sb', transform=transform)

    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    cb1 = plt.colorbar(hdistr, ax=ax, fraction=0.055, aspect=10,  orientation='horizontal', pad=0.1, ticks=lvls[0::2])   
    cb1.set_label(r'$\abs{\dot{\bf{c}}}/\norm{\grad \vb{u}}$')
    cb1.ax.xaxis.set_ticks(lvls, minor=True)
    ax.set_title(titlestr, fontsize=FS, pad=10)

### Plot

geo, prj = sfplt.getprojection(rotation=50+180, inclination=50)

scale = 2.6
fig = plt.figure(figsize=(3/2*1.45*scale,1.00*scale))
gs = fig.add_gridspec(1,3)
gs.update(left=0.04, right=1-0.02, top=0.98, bottom=0.20, wspace=0.4, hspace=0.4)

ax1 = fig.add_subplot(gs[0,0], projection=prj); ax1.set_global(); 
ax2 = fig.add_subplot(gs[0,1], projection=prj); ax2.set_global(); 
ax3 = fig.add_subplot(gs[0,2], projection=prj); ax3.set_global(); 

### Plot 1
ugrad = np.diag([0.5,0.5,-1]) # uniaxial compression along z
plot(ugrad, ax1, titlestr=r'%s Unconfined pure shear'%(r'{\Large \textit{(a)}}\,\,'), speedthres=10e-2)
sfplt.plotcoordaxes(ax1, geo, axislabels='vuxi', color='k')

#### Plot 2
ugrad = np.diag([+1,0,-1]) # confined uniaxial compression along z
plot(ugrad, ax2, titlestr=r'%s Confined pure shear'%(r'{\Large \textit{(b)}}\,\,'))
sfplt.plotcoordaxes(ax2, geo, axislabels='vuxi', color='k')

#### Plot 3
ugrad = np.einsum('i,j',[1,0,0],[0,0,1]) # x-z shear
plot(ugrad, ax3, titlestr=r'%s Simple shear'%(r'{\Large \textit{(c)}}\,\,'))
sfplt.plotcoordaxes(ax3, geo, axislabels='vuxi', color='k')

fout = 'latrot-velocity.png'
print('Saving %s'%(fout))
plt.savefig(fout, transparent=True, dpi=200)

