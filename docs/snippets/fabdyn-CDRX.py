import numpy as np
from specfabpy import specfab as sf
# L=8 truncation is sufficient in this case, but larger L allows a very strong fabric to  
#  develop and minimizes the effect that regularization has on low wavenumber modes (l=2,4)
lm, nlm_len = sf.init(8) 

### Numerics 

Nt = 25   # number of time steps
dt = 0.05 # time step size

### Initial state

nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # state vector
nlm[0,:] = sf.nlm_ideal([0,0,1], 0, L) # normalized single maximum at t=0

### Euler integration
# See Lagrangian parcel demo for more advanced (RK4) integration

for tt in np.arange(1,Nt):

    nlm_prev = nlm[tt-1,:] # previous solution
    Lambda = 1             # CDRX rate factor magnitude
    M = Lambda*sf.M_CDRX(nlm) # CDRX operator 
    nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # Euler step
    nlm[tt,:] = sf.apply_bounds(nlm[tt,:]) # apply spectral bounds if needed
