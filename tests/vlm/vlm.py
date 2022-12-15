# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Given two ODFs of orthogonal axes (coefficieints blm and nlm), calculate ODF of the third axis (coefficients vlm)
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf

res = 40 # latres when plotting ODFs

#----------------------
# Idealized moments for test cases
#----------------------

x,y,z = [1,0,0],[0,1,0],[0,0,1]

### a2

# Single maximum
x2 = np.einsum('i,j',x,x)
y2 = np.einsum('i,j',y,y)
z2 = np.einsum('i,j',z,z)

# Girdle
G2xy = (x2+y2)/2
G2xz = (y2+z2)/2
G2yz = (y2+z2)/2

### a4

# Single maximum
x4 = np.einsum('i,j,k,l',x,x,x,x)
y4 = np.einsum('i,j,k,l',y,y,y,y)
z4 = np.einsum('i,j,k,l',z,z,z,z)

# Girdle
def f_G4pq(p,q):
    pq = np.einsum('i,j',p,q)
    pq2 = pq + pq.T
    p2q2 = np.einsum('i,j',p,p)+np.einsum('i,j',q,q)
    a4 = (np.einsum('i,j,k,l',p,p,p,p) + np.einsum('i,j,k,l',q,q,q,q))/4 \
          + (np.einsum('ij,kl',p2q2,p2q2) + np.einsum('ij,kl',pq2,pq2))/8
    return a4
  
G4xy = f_G4pq(x,y)
G4xz = f_G4pq(x,z)
G4yz = f_G4pq(y,z)

#----------------------

a2_pairs = ((x2,y2), (y2,z2), (G2xy,z2), (G2yz,x2))
a4_pairs = ((x4,y4), (y4,z4), (G4xy,z4), (G4yz,x4))

#----------------------
# Setup figure
#----------------------

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

dpi, scale = 200, 2.5
figsize=(2*scale,2*1.5*scale)
a = 0.00
kwargs_gsupdate = {'left':a, 'right':1-a, 'top':0.99, 'bottom':0.07, 'wspace':0.015*10, 'hspace':0.6}

geo = ccrs.Geodetic()
rot0 = 45
rot = 1.0*rot0 # view angle
inclination = 45 # view angle
prj = ccrs.Orthographic(rot, 90-inclination)

def make_fig():
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(4,3)
    gs.update(**kwargs_gsupdate)
    return fig, gs

def get_axes(gs,ii):
    ax1 = plt.subplot(gs[ii, 0], projection=prj)
    ax2 = plt.subplot(gs[ii, 1], projection=prj)
    ax3 = plt.subplot(gs[ii, 2], projection=prj)
    ax1.set_global(); ax2.set_global(); ax3.set_global()
    return ax1,ax2,ax3

def plot_axes(ax, geo, cax='tab:red'):
    ax.plot([0],[0],  marker=r'$x$', ms=7, c=cax, transform=geo) # x axis
    ax.plot([90],[0], marker=r'$y$', ms=7, c=cax, transform=geo) # y axis
    ax.plot([0],[90], marker=r'$z$', ms=7, c=cax, transform=geo) # z axis

#----------------------
# Test L=2 truncation
#----------------------

print('Calculating for L=2')
lm, nlm_len = sf.init(2)
fig, gs = make_fig()
    
for ii, (a2_b, a2_n) in enumerate(a2_pairs):

    (ax1,ax2,ax3) = get_axes(gs,ii)

    blm = sf.a2_to_nlm(a2_b)
    nlm = sf.a2_to_nlm(a2_n)
    vlm = sf.vlm_L2(blm,nlm)

    plot_ODF(blm, lm, ax=ax1, latres=res, cmap='Greys', cblabel=r'$b/N$ ($L=2$)')
    plot_axes(ax1, geo)

    plot_ODF(nlm, lm, ax=ax2, latres=res, cmap='Greys', cblabel=r'$n/N$ ($L=2$)')
    plot_axes(ax2, geo)

    plot_ODF(vlm, lm, ax=ax3, latres=res, cmap='Greys', cblabel=r'$v/N$ ($L=2$)')
    plot_axes(ax3, geo)

fout = 'vlm-L2.png'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=dpi)

#----------------------
# Test L=4 truncation
#----------------------

print('Calculating for L=4')
lm, nlm_len = sf.init(4)
fig, gs = make_fig()
    
for ii, (a4_b, a4_n) in enumerate(a4_pairs):

    (ax1,ax2,ax3) = get_axes(gs,ii)
    
    blm = sf.a4_to_nlm(a4_b)
    nlm = sf.a4_to_nlm(a4_n)
    vlm = sf.vlm_L4(blm,nlm)

    plot_ODF(blm, lm, ax=ax1, latres=res, cmap='Greys', cblabel=r'$b/N$ ($L=4$)')
    plot_axes(ax1, geo)

    plot_ODF(nlm, lm, ax=ax2, latres=res, cmap='Greys', cblabel=r'$n/N$ ($L=4$)')
    plot_axes(ax2, geo)

    plot_ODF(vlm, lm, ax=ax3, latres=res, cmap='Greys', cblabel=r'$v/N$ ($L=4$)')
    plot_axes(ax3, geo)

fout = 'vlm-L4.png'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=dpi)
    
