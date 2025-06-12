import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(8) 

### Construct an arbitrary fabric to rotate
a2 = np.diag([0, 0, 1]) # arbitrary second-order structure tensor
nlm = np.zeros((nlm_len), dtype=np.complex64) # state vector (array of expansion coefficients)
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients

### Rotate ODF
# Note: assumes L=<12 (rotation for larger L is not implemented)
theta = np.deg2rad(-45) 
phi   = np.deg2rad(45)
nlm_rot1 = sf.rotate_nlm(nlm, theta, 0)    # first rotate around y axis in xz plane
nlm_rot2 = sf.rotate_nlm(nlm_rot1, 0, phi) # next  rotate around z axis in xy plane 
nlm_rot3 = sf.rotate_nlm(nlm_rot2, -theta, -phi) # rotate back
