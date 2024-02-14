# Lagrangian parcel integrator

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation-modes/parcel-trajectory.png#center){: style="width:440px"} 

Specfab includes a high-level integrator for calculating the CPO evolution of a Lagrangian parcel subject to a constant strain-rate (stress) tensor.

```python
import numpy as np
from specfabpy import specfab as sf
from specfabpy import integrator as sfint

L = 12 # expansion series truncation
lm, nlm_len = sf.init(L) 

Nt = 200 # number of integration time steps for below mode of deformation (MOD)
nlm = np.zeros((Nt+1,nlm_len), dtype=np.complex64) # expansion coefficients

"""
For details on the possible mode-of-deformation (MOD) types and options, see
https://nicholasmr.github.io/specfab/deformation-modes
"""

### Pure shear

axis = 2 # axis of compression (0,1,2 = x,y,z)
r = +0.25 # deformation asymmetry in extending directions: r=0 => unconfined, r=+1 or r=-1 => 100% confinement in either of the extending directions
strain_target = -0.95 # deform parcel until this strain target is reached along "axis"
#T = 1/strainrate0 * yr2s # characteristic deformation time-scale of parcel if compressive strain-rate is known (affects only the time vector below)
T = 1

MOD = dict(type='ps', r=r, axis=axis, T=T)

### Simple shear

#plane = 1 # (0,1,2 = yz,xz,xy shear)
#strain_target = np.deg2rad(75) # deform parcel until this target strain (strain angle)
#T = 1 # characteristic deformation time-scale (1/shearrate)

#MOD = dict(type='ss', plane=plane, T=T)

### Integrate

iota   = 1    # plastic spin strength (default: iota=1 for deck-of-cards behaviour)
nu     = 1    # multiplicative factor for default regularization strength (default: nu=1)
Gamma0 = None # DDRX magnitude (None = disabled)
Lambda = None # CDRX magnitude (None = disabled)

nlm[:,:], F, time, ugrad = sfint.lagrangianparcel(sf, MOD, strain_target, Nt=Nt, \
                                iota=iota, nu=nu, Lambda=Lambda, Gamma0=Gamma0)

# nlm (Nt,nlm_len): Harmonic coefficients at each time step
# F (Nt,3,3):       Deformation gradient tensor at each time step
# time (Nt):        Total time at each time step 
# ugrad (3,3):      Velocity gradient

### Auxiliary

strain_ij = np.array([sf.F_to_strain(F[nn,:,:]) for nn in np.arange(Nt)]) # strain_ij tensor


```

