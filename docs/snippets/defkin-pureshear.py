import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) 

axis = 2 # axis of shortening (tauc>0) or lengthening (tauc<0): 0=x, 1=y, 2=z
tauc = 1 # time taken in seconds for parcel to reduce to half (50%) height if tauc>0, or abs(time) taken for parcel to double in height (200%) if tauc<0.
q    = 0 # asymmetry parameter for shortening (if tauc>0) or lengthening (if tauc<0)

tau = tauc/np.log(2)                     # corresponding e-folding time
ugrad = sf.pureshear_ugrad(axis, q, tau) # velocity gradient
D, W = sf.ugrad_to_D_and_W(ugrad)        # strain-rate and spin tensor

t = 1 # some specific time of interest
r   = sf.pureshear_r(tau, t)          # scaling parameter "r" at time t
F   = sf.pureshear_F(axis, q, tau, t) # deformation gradient tensor at time t
eps = sf.F_to_strain(F)               # strain tensor at time t
