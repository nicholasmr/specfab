#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2021

"""
Any ODF(theta,phi) can expanded as a spherical harmonic series, Y_l^m(theta,phi): 
    
    ODF(theta,phi) = sum_{l,m} n_l^m Y_l^m(theta,phi).
    
The second-order structure tensor, a^(2), measures the overlap between the ODF and the lowest-order harmonic modes Y_2^m (for all m=-2,-1,0,1,2):

    a^(2) := <c^2> := <c^2|ODF>  (3x3 matrix)

Given a^(2), the expansion coefficients n_2^m can be reconstructed by writing the ODF in terms of a tensor series expansion, whereby (Siegfried Hess, 2015)

    ODF(theta,phi) = 1/(4*pi) * (1 + 15/3 a^(2):(c^2-I/3) + ...)

Reconstructed the ODF from a^(2) will, however, not represent the true ODF unless all higher-order wavenumber components vanish (c_l^m=0 for l>2, or equivelently a^(i)=0 for i>2).
"""

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

def a2plot(a2):

    # Determine spectral coefs from a^(2)    
    nlm_len = sf.init(2)
    lm      = sf.get_lm(nlm_len)
    clm     = sf.get_a2_to_nlm(nlm_len, a2)

    # Check that re-calculating a2 from clm indeed gives a2 as passed to a2plot()
#    print(1/np.sqrt(4*np.pi), clm)
#    print(a2)
#    print(sf.get_a2_ij(clm))

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

# Plot it
a2plot(np.array(a2))

