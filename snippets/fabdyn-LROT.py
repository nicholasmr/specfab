import numpy as np
from specfabpy import specfab as sf

# L=8 truncation is sufficient in this case, but larger L allows a very strong fabric to  
#  develop and minimizes the effect that regularization has on low wavenumber modes (l=2,4)
lm, nlm_len = sf.init(8) 

### Velocity gradient experienced by parcel

ugrad = np.diag([0.5, 0.5, -1.0]) # uniaxial compression along z-axis
D = (ugrad+np.transpose(ugrad))/2 # symmetric part (strain rate tensor)
W = (ugrad-np.transpose(ugrad))/2 # anti-symmetric part (spin tensor)

### Numerics 

Nt = 25   # number of time steps
dt = 0.05 # time step size

### Initial state

nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # state vector
nlm[0,0] = 1/np.sqrt(4*np.pi) # normalized isotropic state at t=0

### Euler integration

for tt in np.arange(1,Nt):

    nlm_prev = nlm[tt-1,:] # previous solution
    iota, zeta = 1, 0      # "deck of cards" behavior 
    M_LROT = sf.M_LROT(nlm_prev, D, W, iota, zeta) # lattice rotation operator
    M_REG  = sf.M_REG(nlm_prev, D)                 # regularization operator
    M      = M_LROT + M_REG
    nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # Euler step
