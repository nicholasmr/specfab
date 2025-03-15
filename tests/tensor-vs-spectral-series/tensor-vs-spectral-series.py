#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>

"""
Compare tensor and spectral expansion series, based on the same discrete ensemble of c-axes
"""

import copy, sys, os, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp
from scipy.special import factorial as factorial
from fractions import Fraction
from specfabpy import specfab as sf
from specfabpy import plotting as sfplt

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

L = 4
lm, nlm_len = sf.init(L) 

latres=50
lonres=2*50

vcolat = np.linspace(0,   np.pi, latres) # vector
vlon   = np.linspace(0, 2*np.pi, lonres) # vector
lon, colat = np.meshgrid(vlon, vcolat) # gridded (matrix)
lat = np.pi/2-colat
    
#def plot_F(ax, F):
        

def rhat(theta,phi, degree=True):
    if degree:
        theta = np.deg2rad(theta)
        phi   = np.deg2rad(phi)
    return np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])

r = rhat(colat, lon, degree=False)
r2 = np.einsum('aij,bij->abij',r,r)
r4 = np.einsum('aij,bij,cij,dij->abcdij',r,r,r,r)

### Discrete fabric sample

ci = (rhat(30,0),rhat(30,45),rhat(30,90)) #,rhat(30,270))
ci = (rhat(30,0),)

N = len(ci)
a0 = 1 # integral of normalized distribution is 1
a2 = np.zeros((3,3))
a4 = np.zeros((3,3,3,3))
for c in ci:
    a2 += np.einsum('i,j',c,c)/N
    a4 += np.einsum('i,j,k,l',c,c,c,c)/N

if 0: # debug, isotropic
    nlm_iso = np.zeros(nlm_len)
    nlm_iso[0] = 1/np.sqrt(4*np.pi)
    a2 = sf.a2(nlm_iso)
    a4 = sf.a4(nlm_iso)

### Tensorial

a2_irpart = sf.a2_irpart(a2)
a4_irpart = sf.a4_irpart(a2,a4)

qi_num = lambda k: int(factorial(4*k+1))
qi_den = lambda k: int(4**k * factorial(2*k)**2)
qi     = lambda k: 1/(4*np.pi) * qi_num(k)/qi_den(k)
for k2 in (0,1,2):
    print('2k=%i: = %f or '%(k2, qi(k2)), Fraction(qi_num(k2), qi_den(k2)))

F_tensorial = qi(0)*a0 \
            + qi(1)*np.einsum('ab,abij->ij',a2_irpart, r2) \
            + qi(2)*np.einsum('abcd,abcdij->ij',a4_irpart, r4)

#print(F_tensorial, 1/(4*np.pi)) # debug if isotropic fabric 

### Spectral

nlm = sf.a4_to_nlm(a4)
F_spectral = np.real(np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], lon, colat) for ii in np.arange(nlm_len) ], axis=0))

### Plot

geo, prj = sfplt.getprojection(rotation=45, inclination=45)

dpi, scale = 200, 2
figsize=(2.3*scale,1*scale)
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(1,2)
a = 0.1
gs.update(**{'left':a, 'right':1-a, 'top':0.99, 'bottom':0.07, 'wspace':0.02})

ax1 = plt.subplot(gs[0, 0], projection=prj)
ax2 = plt.subplot(gs[0, 1], projection=prj)
ax1.set_global() 
ax2.set_global()

sfplt.plotS2field(ax1, F_spectral, lon, lat) 
sfplt.plotS2field(ax2, F_tensorial, lon, lat) 

#, kwargs_cf={}, \
#            showcb=True, kwargs_cb=kwargs_cb_default, \
#            showgl=True, kwargs_gl=kwargs_gl_default, \
#    ):

plt.savefig('tensor-vs-spectral-series.png', dpi=dpi)

print('err: %.2e'%np.sum(np.abs(F_spectral[:]-F_tensorial[:])))

