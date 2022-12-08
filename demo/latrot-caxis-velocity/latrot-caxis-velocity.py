#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

import copy, sys, code # code.interact(local=locals())
import numpy as np
sys.path.insert(0, '..')
from header import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

### High res

mul = 10
theta = np.linspace(0, np.pi, 20*mul)   
phi   = np.linspace(0, 2*np.pi, 10*mul)
p2, t2 = np.meshgrid(phi, theta)

cv = np.array([np.cos(p2)*np.sin(t2), np.sin(p2)*np.sin(t2), +np.cos(t2)]) 
tv = np.array([np.cos(p2)*np.cos(t2), np.sin(p2)*np.cos(t2), -np.sin(t2)]) 
pv = np.array([-np.sin(p2), np.cos(p2), 0*p2]) 

### Low res

mul = 1.0
theta_ = np.linspace(0, np.pi, int(10*mul))   
phi_   = np.linspace(0, 2*np.pi, int(20*mul))
p2_, t2_ = np.meshgrid(phi_, theta_)

cv_ = np.array([np.cos(p2_)*np.sin(t2_), np.sin(p2_)*np.sin(t2_), +np.cos(t2_)]) 
tv_ = np.array([np.cos(p2_)*np.cos(t2_), np.sin(p2_)*np.cos(t2_), -np.sin(t2_)]) 
pv_ = np.array([-np.sin(p2_), np.cos(p2_), 0*p2_]) 

### Functions

def plot(ugrad, ax, titlestr='', speedthres=0):

    ### Contour high res

    D, W = (ugrad + ugrad.T)/2, (ugrad - ugrad.T)/2

    cdot = np.einsum('ij,jkl->ikl', W, cv) - (np.einsum('ij,jkl->ikl', D, cv) - np.einsum('ijk,ljk,lm,mjk->ijk', cv,cv,D,cv))
         
    u = np.einsum('ijk,ijk->jk', cdot, pv)
    v = np.einsum('ijk,ijk->jk', cdot, tv)
    speed = np.linalg.norm(cdot, axis=0)
    u[speed<speedthres] = np.nan

    lvls = np.linspace(0,0.8,5)
    x, y = np.rad2deg(p2), np.rad2deg(t2 - np.pi/2) 
    F = speed/np.linalg.norm(ugrad)
    hdistr = ax.contourf(x,y, F, transform=ccrs.PlateCarree(), extend='max', cmap='YlOrBr', levels=lvls)
#    ax.streamplot(x,y,u,v, color='0', linewidth=0.9, density=[0.45,0.65], arrowsize=1,  transform=geo)

    ### Quiver, low res

    cdot_ = np.einsum('ij,jkl->ikl', W, cv_) - (np.einsum('ij,jkl->ikl', D, cv_) - np.einsum('ijk,ljk,lm,mjk->ijk', cv_,cv_,D,cv_))
         
    u_ = np.einsum('ijk,ijk->jk', cdot_, pv_)
    v_ = np.einsum('ijk,ijk->jk', cdot_, tv_)
    speed_ = np.linalg.norm(cdot_, axis=0)
    u_[speed_<speedthres] = np.nan

    cquiv = 'k' # #b2182b       
#    cquiv = '#a50f15' # #b2182b
    lvls = np.linspace(0,0.8,5)
    x_, y_ = np.rad2deg(p2_), np.rad2deg(t2_ - np.pi/2) 
    QV1 = ax.quiver(x_,y_,u_,v_, scale=7.5, width=0.012, color=cquiv, transform=geo)
    plt.quiverkey(QV1, 0.95, 0.00, 1, r'$\dot{\vb{c}}$', labelpos='E', coordinates='axes', labelsep=0.05)

    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    cb1 = plt.colorbar(hdistr, ax=ax, fraction=0.055, aspect=10,  orientation='horizontal', pad=0.1, ticks=lvls[0::2])   
    cb1.set_label(r'$\abs{\dot{\bf{c}}}/\norm{\grad \vb{u}}$')
    cb1.ax.xaxis.set_ticks(lvls, minor=True)
    ax.set_title(titlestr, fontsize=FS, pad=10)

### Plot

inclination = 45 # view angle
rot0 = -90
rot = -45 + rot0 # view angle

prj = ccrs.Orthographic(rot, 90-inclination)
geo = ccrs.RotatedPole()

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

### Plot 1
ugrad = np.diag([0.5,0.5,-1]) # Unconfined UC in z
plot(ugrad, ax1, titlestr=r'%s Unconfined pure shear'%(r'{\Large \bf a}\,\,'), speedthres=10e-2)

#### Plot 2
ugrad = np.diag([+1,0,-1]) 
plot(ugrad, ax2, titlestr=r'%s Confined pure shear'%(r'{\Large \bf b}\,\,'))

#### Plot 3
ugrad = np.array([[0,0,1], [0,0,0], [0,0,0]])
plot(ugrad, ax3, titlestr=r'%s Simple shear'%(r'{\Large \bf c}\,\,'))

#---------

for ax in axlist:
    ax.plot([0],[90], c='0.3', marker=r'$z$', ms=8, transform=geo)
    ax.plot([0],[0], c='0.3', marker=r'$x$', ms=8, transform=geo)
    ax.plot([90],[0], c='0.3', marker=r'$y$', ms=8, transform=geo)

#---------

fout = 'latrot-caxis-velocity.png'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=200)

