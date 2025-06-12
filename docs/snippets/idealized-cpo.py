import numpy as np
from specfabpy import specfab as sf
L = 8
lm, nlm_len = sf.init(L) 

m = [0,0,1] # symmetry axis of orientation distribution 
# colat=0    => unidirectional distribution
# colat=pi/2 => planar distribution
# ...anything in between is a small circle distribution
colat = 0 
nlm = sf.nlm_ideal(m, colat, L) # note: only l<=12 coefs are determined even if L>12
