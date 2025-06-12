import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4)
nlm = np.zeros((nlm_len), dtype=np.complex64) # state vector (array of expansion coefficients)

### Physical parameters (note eps=eps0*epsr, mu=mu0*mur)
epsr_c = 3.17        # relative permittivity of a single grain parallel to symmetry axis (c) 
epsr_a = 3.17-0.034  # relative permittivity of a single grain perpendicular to symmetry axis (a)
mur = 1              # relative permeability of a single grain

### CPO from second-order structure tensor
p = np.array([0,0,1]) # preferred c-axis direction
a2 = np.einsum('i,j', p,p) # = deltafunc(r-p) 
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients

### Propagation directions of interest
theta, phi = np.deg2rad([0,90,]), np.deg2rad([0,0,]) # theta=colatitude, phi=longitude

### Calculate phase velocities
Vi = sf.Vi_electromagnetic_tranisotropic(nlm, epsr_c, epsr_a, mur, theta,phi) # fast and slow phase velocities are V_S1=vi[0,:], V_S2=vi[1,:]
