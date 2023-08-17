# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

"""
Header file for common routines etc.
"""

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
from scipy.spatial.transform import Rotation
import quaternion as qt # pip3 install numpy-quatern (if fails try pip3 install numpy==1.23.1)
import pandas as pd

import matplotlib.pyplot as plt

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt

norm = 1/np.sqrt(4*np.pi)

# fwd/inv transform pair for n_2^0 <--> a_zz^(2)
def n20_to_azz(x): return 1/3*(1+2/np.sqrt(5)*x)  # fwd
def azz_to_n20(x): return (x - 1/3)/(2/np.sqrt(5)) # inv

def get_m_angles(a2, Ilam1=-1):

    """
    Get a^(2) eigenvector angles need to rotate (measured) ODFs into vertical frame
        If Ilami=-1 => rotation angle for the largest eigenvalue direction (used for single max)
        If Ilami=0  => rotation angle for the smallest eigenvalue direction (used for girdles)
    """
    
    lami, mi = np.linalg.eig(a2)
    I = np.argsort(lami)[Ilam1]
    m1 = mi[:,I] 
    _, m1_colat, m1_lon = sfdsc.cart2sph(m1)
    return (m1_colat, m1_lon, m1) 
    
def plot_points(ax, latd, lond, **kwargs):

    colat, lon = sfdsc.lat2colat(np.deg2rad(latd)), np.deg2rad(lond)
    v = sfdsc.sph2cart(colat, lon)
    sfplt.plotS2point(ax, +v, ls='none', **kwargs)
    sfplt.plotS2point(ax, -v, ls='none', **kwargs) # antipodal symmetric

def plot_trajectory(ax, nlm, arrpos=None, label=None, ls='-', lw=1.75, c='k', hwmul=1, zorder=14,  endmarker=False, mse=7):

    x, y = np.real(nlm[:,sf.I20]), np.real(nlm[:,sf.I40])
    h, = ax.plot(x,y, ls=ls,lw=lw,c=c, label=label, zorder=zorder)
   
    if arrpos is not None:
        jj, djj = arrpos, 4
        xa0,ya0, xa1,ya1 = x[jj],y[jj], x[jj+djj],y[jj+djj]
        dx, dy = xa1-xa0, ya1-ya0
        arrlen = np.sqrt(dx**2 + dy**2)
        dx *= 0.5e-2/arrlen
        dy *= 0.5e-2/arrlen
        ax.arrow(xa0, ya0, dx, dy, shape='full', lw=1, color=c,  head_width=.02*2.2*hwmul/normfac, head_length=0.035*0.8/normfac, zorder=zorder)
    
    if endmarker: ax.plot(x[-1],y[-1], marker='o', ms=mse, ls='none', c=c, label=None, zorder=zorder) # fillstyle='none',
    
    return h

def load_sample(fname, expr):

    fullfname = 'data/%s/%s'%(expr['path'],fname)
    print('\n*** Loading %s'%(fullfname))

    df = pd.read_csv(fullfname, skiprows=expr['skiprows'], sep=expr['sep']) # header=expr['headerlen'], 

    # For SPICE files remove rows with "unrecognized grains"
    m = ~df.astype(str).apply(lambda x: x.str.contains('grain', regex=False)).any(axis=1)
    df = df[m]            

    # Quaternion coordinates?        
    if expr['coords'] == 'quaternion':
        # Typical for processed EBSD data in MTEX (MATLAB)
        qs_comps = df.to_numpy()[:,:4] 
        print('- using quaternion representation :: shape(qs_comps) =', np.shape(qs_comps))
        qs = qt.as_quat_array(qs_comps)
        sphcoords = qt.as_spherical_coords(qs)
        qcolat, qlon = sphcoords[:,0], sphcoords[:,1] # colat [0;pi], lon [0;2pi]
        qlat = sfdsc.colat2lat(qcolat)

    # Cartesian coordinate vector?            
    elif expr['coords'] == 'vector': 
        vj = df.to_numpy()
        print('- using vector representation :: shape(vj) =', np.shape(vj))
        qlat, qcolat, qlon = sfdsc.cart2sph(vj)

    # Spherical coordinates?        
    else: 
        sphcoords = df.to_numpy()
        print('- using spherical coordinate representation :: shape(sphcoords) =',np.shape(sphcoords))
        qlat = np.deg2rad(sphcoords[:,expr['I_lat']].astype(np.float64))
        qlon = np.deg2rad(sphcoords[:,expr['I_azi']].astype(np.float64))
        qcolat = sfdsc.lat2colat(qlat)
        if expr['iscolat']: qlat, qcolat = qcolat, qlat # swap

    ### Determine ODF from discrete c-axis ensemble
    
    lm, nlm_len = sf.init(6)
    caxes = sfdsc.sph2cart(qcolat, qlon)
    nlm = sfdsc.vi2nlm(caxes, L=6) 

    ### Rotated frame 
       
    # Preferred direction (postulated approx. symmetry axis)
    Ilam1 = 0 if expr['type'] == 'ue' else 2 # sym. axis = largest eig dir. for single max, smallest eig dir. for girdle
    (m_colat, m_lon, _) = get_m_angles(sf.a2(nlm), Ilam1=Ilam1)
    #m_colat = 0 # debug: only horiz rotation

    # True nlm array from which spectral coefs are derived for correlation        
    print('- old m vector colat, lon = %.1f, %.1f (deg.) '%(np.rad2deg(m_colat), np.rad2deg(m_lon)))
    nlmr = sf.rotate_nlm(nlm, 0, -m_lon)    # rotation is the composition of two rotations 
    nlmr = sf.rotate_nlm(nlmr, -m_colat, 0) # ... second rotation

    # verify that rotated a2 frame has m vertical
    (mnew_colat, mnew_lon, mnew) = get_m_angles(sf.a2(nlmr))
    print('- new m vector colat, lon = %.1f, %.1f (deg.) '%(np.rad2deg(mnew_colat), np.rad2deg(mnew_lon)))
    print('- new m vector = (%.2e, %.2e, %.2e)'%(mnew[0],mnew[1],mnew[2]))
    
    # Determine c-axes lat+lon in rotated frame
    Ry = Rotation.from_euler('y', -m_colat).as_matrix()
    Rz = Rotation.from_euler('z', -m_lon).as_matrix()
    caxes = np.array([np.einsum('ij,jk,k',Ry,Rz,c) for c in caxes])
    qlatr, qcolatr, qlonr = sfdsc.cart2sph(caxes)

    ### Return 
    
    q  = (qlat, qcolat, qlon)
    qr = (qlatr, qcolatr, qlonr)
    m  = (sfdsc.colat2lat(m_colat), m_colat, m_lon)
    mew = (sfdsc.colat2lat(mnew_colat), mnew_colat, mnew_lon)
    return (q, qr, m, mnew, caxes, nlm, nlmr, lm)
    
