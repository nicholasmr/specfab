import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) # L=4 is sufficient here
nlm = np.zeros((nlm_len), dtype=np.complex64) 
blm = np.zeros((nlm_len), dtype=np.complex64) 

### CPO (nlm,blm) from fourth-order structure tensors
vn = np.array([0,0,1]) # preferred slip-normal direction
vb = np.array([1,0,0]) # preferred slip direction
a4_n = np.einsum('i,j,k,l', vn,vn,vn,vn) # n(r) = deltafunc(r-vn) 
a4_b = np.einsum('i,j,k,l', vb,vb,vb,vb) # b(r) = deltafunc(r-vb) 
nlm[:sf.L4len] = sf.a4_to_nlm(a4_n) # l<=4 expansion coefficients
blm[:sf.L4len] = sf.a4_to_nlm(a4_b) # l<=4 expansion coefficients

### Physical parameters (SI units)
rho = 3355 # density of olivine
alpha = 1 # Voigt--Reuss weight; only alpha=1 supported for now
Cij = (320.5e9, 196.5e9, 233.5e9,  64.0e9, 77.0e9, 78.7e9,  76.8e9, 71.6e9, 68.15e9) # Abramson (1997) parameters (C11,C22,C33,C44,C55,C66,C23,C13,C12)
Lame_grain = sf.Cij_to_Lame_orthotropic(Cij) # Lame parameters (lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3)
# Note that the above ordering of Lame/Cij parameters assume an A-type fabric; that is, (blm,nlm,vlm) refer to the distibutions of (m1',m2',m3') axes, respectively.
# If concerned with another fabric type, the components can easily be re-ordered:
#Lame_grain = sf.Lame_olivine_A2X(Lame_grain, 'B') # B-type Lame paremeters
#Lame_grain = sf.Lame_olivine_A2X(Lame_grain, 'C') # C-type Lame paremeters

### Propagation directions of interest
theta, phi = np.deg2rad([90,70,]), np.deg2rad([0,10,]) # theta=colatitude, phi=longitude

### Calculate phase velocities
vlm = 0*nlm # estimate vlm from (blm,nlm) by passing zero array
Vi = sf.Vi_elastic_orthotropic(blm,nlm,vlm, alpha, Lame_grain, rho, theta,phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
