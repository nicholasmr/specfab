import numpy as np
from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
lm, nlm_len = sf.init(8) 

"""
Generate sample of N grains with random c-axis orientations and sizes
"""

N     = 1000 # number of grains
lon   = np.random.uniform(0, 2*np.pi, N)       # uniform longitude
colat = np.arccos(np.random.uniform(-1, 1, N)) # uniform colatitude
caxes = sfdsc.sph2cart(colat, lon)
weights = np.random.uniform(0.5, 1, N) # grain weights (typically grain sizes)
weights /= np.linalg.norm(weights)     # normalize weights

"""
Determine state vector nlm (array of harmonic expansion coefficients)
"""

### Do it all yourself

a6 = np.zeros((3,3,3,3,3,3))
for c in caxes:
    a6 += np.einsum('i,j,k,l,m,n', c,c,c,c,c,c) # sixth outer product of c-axis with itself
a6 /= N # normalize by number of grains if grain sizes are not known
nlm = np.zeros((nlm_len), dtype=np.complex64) # state vector (expansion coefficients)
nlm[:sf.L6len] = sf.a6_to_nlm(a6) # l<=6 expansion coefficients
print(nlm)

### Alternatively, use built-in function ri_to_nlm(caxes, weights, L)

weights = np.ones(N)/N # equal wight to each grain, should sum to 1
nlm[:sf.L6len] = sf.ri_to_nlm(caxes, weights, 6) # l<=6 expansion coefficients
print(nlm)
