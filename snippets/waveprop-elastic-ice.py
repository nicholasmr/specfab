import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) # L=4 is sufficient here
nlm = np.zeros((nlm_len), dtype=np.complex64) # state vector (array of expansion coefficients)

### CPO from fourth-order structure tensor
p = np.array([0,0,1]) # preferred c-axis direction
a4 = np.einsum('i,j,k,l', p,p,p,p) # = deltafunc(r-p) 
nlm[:sf.L4len] = sf.a4_to_nlm(a4) # l<=4 expansion coefficients

### Physical parameters (SI units)
rho = 917 # density of ice
Cij = (14.060e9, 15.240e9, 3.060e9, 7.150e9, 5.880e9) # Bennett (1968) parameters (C11,C33,C55,C12,C13)
Lame_grain = sf.Cij_to_Lame_tranisotropic(Cij) # Lame parameters (lam,mu,Elam,Emu,Egam)
alpha = 0.5 # Voigt--Reuss weight, where 0.5 = Hill average

### Propagation directions of interest
theta, phi = np.deg2rad([90,70,]), np.deg2rad([0,10,]) # theta=colatitude, phi=longitude

### Calculate phase velocities
Vi = sf.Vi_elastic_tranisotropic(nlm, alpha, Lame_grain, rho, theta,phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
