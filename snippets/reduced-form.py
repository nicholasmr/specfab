import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(6)

"""
Construct an arbitrary fabric state
"""
a2 = np.diag([0.1,0.2,0.7]) # some second-order structure tensor
nlm = np.zeros((nlm_len), dtype=np.complex64) # full form of state vector
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients from a2
print('original:', nlm)

"""
Get reduced form of state vector, rnlm
"""
rnlm_len = sf.get_rnlm_len() 
rnlm = np.zeros((rnlm_len), dtype=np.complex64) # reduced state vector
rnlm[:] = sf.nlm_to_rnlm(nlm, rnlm_len) # reduced form
print('reduced:', rnlm)

"""
Recover full form (nlm) from reduced form (rnlm)
"""
nlm[:] = sf.rnlm_to_nlm(rnlm, nlm_len)
print('recovered:', nlm)
