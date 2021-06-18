#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2021

import sys, os, code, copy # code.interact(local=locals())
import numpy as np
import scipy.special as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

from specfabpy import specfabpy as sf

#------------------
# Plot resolution and options
#------------------

inclination = 60 # View angle
rot = -90 # View angle
prj, geo = ccrs.Orthographic(rot, 90-inclination), ccrs.Geodetic()

latres = 20
lonres = 2*latres
phi   = np.linspace(0, 2*np.pi, lonres) # LON
theta = np.linspace(0,   np.pi, latres) # CO-LAT 
phi, theta = np.meshgrid(phi, theta)            
lon, colat = phi, theta
lat = np.pi/2-colat

lvls = np.linspace(0, 0.4, 6) # Contour levels

#------------------
# Plot function
#------------------

def a4plot(a2,a4):

    # Determine spectral coefs from a^(2) and a^(4)
    nlm_len = sf.init(4)
    lm      = sf.get_lm(nlm_len)
    clm     = sf.get_a4_to_nlm(nlm_len, a2, a4)

    # Check that re-calculating a2,a4 from clm indeed gives arguments passed to a4plot()
    #print(1/np.sqrt(4*np.pi), clm)
    #print(a2)
    #print(sf.get_a2_ij(clm))
    #print('---')
    #print('max[abs[delta a^(4)]] =', np.amax(np.abs(a4-sf.get_a4_ijkl(clm))))

    # Discretize ODF 
    ODF = np.sum([ clm[ii]*sp.sph_harm(m, l, phi,theta) for ii,(l,m) in enumerate(lm.T) ], axis=0)
    ODF = np.real(ODF) # Any imag values are rounding errors.
    
    # Plot ODF
    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0, 0], projection=prj)
    ax.set_global();
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    h = ax.contourf(np.rad2deg(lon),np.rad2deg(lat), ODF, transform=ccrs.PlateCarree(), levels=lvls, extend='both', cmap=cmap)
    h.cmap.set_under('#fee0d2')
    cb = plt.colorbar(h, ax=ax)   
    g = ax.gridlines(crs=ccrs.PlateCarree())
    cb.set_label(r'$\mathrm{ODF}\left( a^{(2)} \right)$')
    plt.show()

#------------------
# Examples
#------------------

## a^(2)

# Single maximum
k1 = 0.99 # Largest eigen value: 1 >= k_1 >= 1/3
k23 = (1-k1)/2
#a2 = [[k1,0,0],[0,k23,0],[0,0,k23]] # max in x
#a2 = [[k23,0,0],[0,k1,0],[0,0,k23]] # max in y
a2 = [[k23,0,0],[0,k23,0],[0,0,k1]] # max in z

# Girdle
k12 = 0.5 # Largest eigen values: 0.5 >= k_1=k_2 >= 1/3
k3 = 1-2*k12
#a2 = [[k12,0,0],[0,k12,0],[0,0,k3]] # x-y girdle
#a2 = [[k12,0,0],[0,k3,0],[0,0,k12]] # x-z girdle

a2 = np.array(a2)

## a^(4) 

#   Note that a^(4) requires subtraces=1 *and* contains contributions from a^(2) (in the sense that a^(4) depends on n_l^m for l<=4, and a^(2) on n_l^m for l<=2).
#   Hence, the entries cannot be freely set for our synthetic example in analogy to a^(2) above. 
#   Insted, we specify the independent contributions to a^(4) by setting the corresponding n_4^m coefficients, which easily allows for a set of "self-consistent" tensors a^(2) and a^(4).

nlm_len = sf.init(4)
nlm = sf.get_a2_to_nlm(nlm_len, a2) # get n_2^m spectral coefs from a^(2)
nlm[11] = 1 # Set some higher-order anisotropy (indices 6 <= i <= 14 are the n_4^m coefs)
a4  = sf.get_a4_ijkl(nlm) 

## Plot it
a4plot(a2, a4)

