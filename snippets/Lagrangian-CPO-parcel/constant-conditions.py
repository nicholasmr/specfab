import numpy as np
from specfabpy import specfab as sf
from specfabpy import integrator as sfint

L = 12 # expansion series truncation
lm, nlm_len = sf.init(L) 

Nt = 200 # number of integration time steps taken
nlm = np.zeros((Nt+1,nlm_len), dtype=np.complex64) # expansion coefficients

### Pure shear

axis = 2     # axis of compression (0,1,2 = x,y,z)
tau  = 1     # characteristic e-folding deformation time
q    = +0.25 # asymmetry in directions of extension: q=0 => unconfined, q=+1 or q=-1 => 100% confinement in either of the extending directions
strain_target = -0.95 # deform parcel until this strain target is reached along "axis"

DK = dict(type='ps', q=q, axis=axis, tau=tau) # deformation kinematics

### Simple shear

#plane = 1 # (0,1,2 = yz,xz,xy shear)
#strain_target = np.deg2rad(75) # deform parcel until this target strain (strain angle)
#tau = 1 # characteristic deformation time-scale (1/shearrate)

#DK = dict(type='ss', plane=plane, tau=tau)

### Integrate

kw_dyn = dict(
    iota   = 1,    # plastic spin strength (default: iota=1 for deck-of-cards behaviour)
    nu     = 1,    # multiplicative factor for default regularization strength (default: nu=1)
    Gamma0 = None, # DDRX magnitude (None = disabled), assumes stress is coaxial to strain rate
    Lambda = None, # CDRX magnitude (None = disabled)
)

nlm[:,:], Fi, ti, ugrad = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, **kw_dyn)

"""
Outputs are:
------------
nlm (Nt,nlm_len): Harmonic coefficients at each time step
Fi (Nt,3,3):      Deformation gradient tensor at each time step
ti (Nt):          Total time at each time step 
ugrad (3,3):      Velocity gradient
"""

### Auxiliary

strain = np.array([sf.F_to_strain(Fi[nn]) for nn in np.arange(Nt)]) # strain tensor
