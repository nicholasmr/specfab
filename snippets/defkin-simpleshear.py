import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) 

plane = 1 # plane of shear: 0=yz, 1=xz, 2=xy
tau   = 1 # time taken in seconds for parcel to a reach shear strain of 1 (45 deg shear angle)

ugrad = sf.simpleshear_ugrad(plane, tau) # velocity gradient
D, W = sf.ugrad_to_D_and_W(ugrad)        # strain-rate and spin tensor

t = 1 # some specific time of interest
gamma = sf.simpleshear_gamma(tau, t)    # shear angle at time t
F     = sf.simpleshear_F(plane, tau, t) # deformation gradient tensor at time t
eps   = sf.F_to_strain(F)               # strain tensor at time t
