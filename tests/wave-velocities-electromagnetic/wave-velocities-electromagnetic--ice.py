# Nicholas Rathmann <rathmann@nbi.ku.dk>

"""
Test wave velocities predicted for a single ice crystal
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import discrete as sfdsc

L = 4
lm, nlm_len = sf.init(L)

### CPO with perfectly aligned grains (same behaviour as a single grain)

x,y,z = np.eye(3)

r = [0,0,1]
nlm = sf.nlm_ideal(r, 0, L) # [100] axis aligned with x
print(sf.a2(nlm))

#nlm = sf.a2_to_nlm(np.diag([0.0,0.3,0.7]))

### Materials and wave params

epsa = 3.17-0.034*10
epsc = 3.17

#epsa, epsc = 3.2, 3.2

mu = 1

### Calculate wave speeds

def get_sphr_coords(v):
    x,y,z, r = v[0],v[1],v[2], 1
    theta = np.arccos(z/r) # co-latitude
    phi = np.arctan2(y, x) 
    return theta, phi 
        
khat_names = ['x','y','z']
khat = (x,y,z)
        
angdir = np.array([ np.array(sfdsc.cart2sph(kh))[[1,2]] for kh in khat])
Vi = sf.Vi_electromagnetic_tranisotropic(nlm, epsc, epsa, mu, angdir[:,0],angdir[:,1])

Vi *= 1e-6 # m/microsecond
print('*** units are m/microsecond ***')

for ii in range(len(khat)):
    print('----')
    print('%s => (colat, lon) = (%.0f, %.0f) deg.'%(khat_names[ii], np.rad2deg(angdir[ii,0]), np.rad2deg(angdir[ii,1])))
    print('(vS1, vS2) = (%.0f, %.0f)'% (Vi[0,ii],Vi[1,ii]) )
    
#    print('(vS1, vS2) = (%.0f, %.0f, %.0f, %.0f)'% (Vi[0,ii],Vi[1,ii], Vi[2,ii],Vi[3,ii]) )

