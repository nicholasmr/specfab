import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(8) 

### Replace with your own array/list of measured c-axes
# caxes = [[c1x,c1y,c1z], [c2x,c2y,c2z], ...] 

### Determine sixth-order structure tensor
a6 = np.zeros((3,3,3,3,3,3))
for c in caxes
    a6 += np.einsum('i,j,k,l,m,n', c,c,c,c,c,c) # sixth outer product of c-axis with itself
a6 /= len(caxes) # normalize by number of c-axes (grains) if grain sizes are not known

### Determine state vector
nlm = np.zeros((nlm_len), dtype=np.complex64) # state vector (array of expansion coefficients)
nlm[:sf.L6len] = sf.a6_to_nlm(a6) # l<=6 expansion coefficients
