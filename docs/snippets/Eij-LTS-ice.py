import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(8) 

### Synthetic unidirectional CPO (all c-axes aligned in z-direction)
m = np.array([0,0,1]) 
a4 = np.einsum('i,j,k,l', m,m,m,m) # fourth-order structure tensor
nlm = np.zeros((nlm_len), dtype=np.complex64)
nlm[:sf.L4len] = sf.a4_to_nlm(a4) # corresponding expansion coefficients

### Enhancement-factor basis
(e1,e2,e3, eigvals) = sf.frame(nlm, 'e') # enhancement factors are w.r.t. a^(2) basis (i.e. eigenenhancements)
#(e1,e2,e3) = np.eye(3) # enhancement factors are w.r.t. Cartesian basis (x,y,z)

### Transversely isotropic monocrystal parameters for ice (Rathmann & Lilien, 2021)
n_grain   = 1        # power-law exponent: n=1 => linear grain rheology, nonlinear (n>1) is unsupported.
Eij_grain = (1, 1e3) # grain eigenenhancements (Ecc,Eca) for compression along c-axis (Ecc) and for shear parallel to basal plane (Eca)
alpha     = 0.0125   # Taylor--Sachs weight

### Calculate enhancement factors w.r.t. (e1,e2,e3)
Eij = sf.Eij_tranisotropic(nlm, e1,e2,e3, Eij_grain,alpha,n_grain) # Eij=(E11,E22,E33,E23,E13,E12)
