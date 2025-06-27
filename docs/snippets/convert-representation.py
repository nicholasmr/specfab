import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(8)
nlm = np.zeros((nlm_len), dtype=np.complex64) # state vector (expansion coefficients)

"""
a2 to nlm
"""

a2 = np.diag([0.0,0.25,0.75]) # some second-order structure tensor
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients
a2 = sf.a2(nlm) # nlm back to a2
print('a2 is: ', a2)

"""
a4 to nlm
"""

p = np.array([0,0,1]) # unidirectional CPO
a4 = np.einsum('i,j,k,l', p,p,p,p) # = deltafunction(r-p) 
nlm[:sf.L4len] = sf.a4_to_nlm(a4) # l<=4 expansion coefficients
a4 = sf.a4(nlm) # nlm back to a4 
print('a4 is: ', a4)

"""
a6 to nlm
"""

p = np.array([0,0,1]) # unidirectional CPO
a6 = np.einsum('i,j,k,l,m,n', p,p,p,p,p,p) # = deltafunction(r-p) 
nlm[:sf.L6len] = sf.a6_to_nlm(a6) # l<=6 expansion coefficients
a6 = sf.a6(nlm) # nlm back to a6
print('a6 is: ', a6)
