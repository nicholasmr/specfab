# Nicholas Rathmann <rathmann@nbi.ku.dk>

"""
Test wave velocities predicted for a single olivine crystal
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import discrete as sfdsc

L = 4
lm, nlm_len = sf.init(L)

### CPO with Perfectly aligned grains (same behaviour as a single grain)

x, y, z = np.eye(3)
nlm_r1 = sf.nlm_ideal(x, 0, L) # [100] axis aligned with x
nlm_r2 = sf.nlm_ideal(y, 0, L) # [010] axis aligned with y
nlm_r3 = sf.nlm_ideal(z, 0, L) # [001] axis aligned with z

### Olivine monocrystal elastic parameters

rho = sfconst.olivine['density']
Cij = sfconst.olivine['elastic']['Jacobsen2008'] # (C11,C22,C33,C44,C55,C66,C23,C13,C12)
Lame_grain = sf.Cij_to_Lame_orthotropic(Cij) # Lame parameters (lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3)
alpha = 1 # Voigt--Reuss weight; only alpha=1 supported for now

### Calculate phase velocities

def get_sphr_coords(v):
    x,y,z, r = v[0],v[1],v[2], 1
    theta = np.arccos(z/r) # co-latitude
    phi = np.arctan2(y, x) 
    return theta, phi 
        
khat_names = ['x (a)','y (b)','z (c)']
khat = (x,y,z)
        
angdir = np.array([ np.array(sfdsc.cart2sph(kh))[[1,2]] for kh in khat])
Vi = sf.Vi_elastic_orthotropic(nlm_r1,nlm_r2,nlm_r3, alpha, Lame_grain, rho, angdir[:,0],angdir[:,1])

for ii in range(len(khat)):
    print('----')
    print('%s => (colat, lon) = (%.0f, %.0f) deg.'%(khat_names[ii], np.rad2deg(angdir[ii,0]), np.rad2deg(angdir[ii,1])))
    print('(vS1, vS2, vP) = (%.0f, %.0f, %.0f) \nvS1-vS2 = %1.f'% (Vi[0,ii],Vi[1,ii],Vi[2,ii], Vi[0,ii]-Vi[1,ii]) )

