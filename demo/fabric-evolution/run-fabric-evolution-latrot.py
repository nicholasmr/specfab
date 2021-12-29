# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020-2021

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp
from netCDF4 import Dataset

sys.path.insert(0, '..')
from specfabpy import specfabpy as sf # To use specfabpy compile the specfab Python module by running "make specfabpy"

#----------------------
# Input arguments
#----------------------

if len(sys.argv) != 2:
    print('usage: python3 %s ugrad'%(sys.argv[0]))
    sys.exit(0)
    
arg_exp = sys.argv[1] # Velocity gradient experiment

#----------------------
# Numerics
#----------------------

Nt = 50 # Number of time steps
dt = 0.0782404601085629 # Time-step size (gives a vertical strain of -0.98 for experiment "uc_zz")
L = 6 # Spectral truncation (4<=L<=8)

#----------------------
# Grain parameters
#----------------------

Eca_lin,  Ecc_lin,  alpha_lin  = sf.Eca_opt_lin,  sf.Ecc_opt_lin,  sf.alpha_opt_lin  # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
Eca_nlin, Ecc_nlin, alpha_nlin = sf.Eca_opt_nlin, sf.Ecc_opt_nlin, sf.alpha_opt_nlin # Optimal n'=3 (nlin) grain parameters (Rathmann et. al, 2021)

#----------------------
# Velocity gradient
#----------------------

# "ugrad" is the specified (constant) velocity gradient.

ugrad = np.diag([0.5, 0.5, -1.0]) 

if arg_exp == 'uc_xx': ugrad = np.diag([-1, .5, .5])
if arg_exp == 'uc_yy': ugrad = np.diag([.5, -1, .5])
if arg_exp == 'uc_zz': ugrad = np.diag([.5, .5, -1])

if arg_exp == 'ue_xx': ugrad = -1*np.diag([-1, .5, .5])
if arg_exp == 'ue_yy': ugrad = -1*np.diag([.5, -1, .5])
if arg_exp == 'ue_zz': ugrad = -1*np.diag([.5, .5, -1])

if arg_exp == 'cc_zx': ugrad = np.diag([1, 0, -1])
if arg_exp == 'cc_zy': ugrad = np.diag([0, 1, -1])
if arg_exp == 'cc_yx': ugrad = np.diag([1, -1, 0])

if arg_exp == 'ss_xz': ugrad = np.array([[0,0,1], [0,0,0], [0,0,0]])
if arg_exp == 'ss_xy': ugrad = np.array([[0,1,0], [0,0,0], [0,0,0]])
if arg_exp == 'ss_yz': ugrad = np.array([[0,0,0], [0,0,1], [0,0,0]])

if arg_exp == 'rr_xz': ugrad = np.array([[0,0,+1], [0,0,0], [-1,0,0]])
if arg_exp == 'rr_xy': ugrad = np.array([[0,+1,0], [-1,0,0], [0,0,0]])
if arg_exp == 'rr_yz': ugrad = np.array([[0,0,0], [0,0,+1], [0,-1,0]])

eps = (ugrad+np.transpose(ugrad))/2 # Symmetric part (strain-rate)
omg = (ugrad-np.transpose(ugrad))/2 # Anti-symmetric part (spin)

#----------------------
# Initialize model
#----------------------

# The (l,m) values corresponding to the coefficients in "nlm".
# nlm_len is the number of fabric expansion coefficients (degrees of freedom).
lm, nlm_len = sf.init(L) 
nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # The expansion coefficients

# Initial fabric state
if arg_exp[0:2] == 'rr':
    # Init with some anisotropy
    nlm[0,0] = np.sqrt(1/2 + 0j)
    nlm[0,2] = np.sqrt(1/2 + 0j)
else:
    # Init with isotropy
    nlm[:,0] = 1/np.sqrt(4*np.pi) # Normalized such that N(t=0) = 1

#----------------------
# Integrate
#----------------------

if arg_exp == 'uc_zz' or arg_exp == 'cc_zx' or arg_exp == 'ue_zz': 
    te = 1/eps[-1,-1]
    print(r'Integrating until strain_zz = %.3f'%(np.exp(Nt*dt/te)-1))

print('Numerics: Nt=%i, dt=%f, L=%i (nlm_len=%i)'%(Nt,dt,L,nlm_len))

# Euler integration scheme
for tt in np.arange(1,Nt):

    nlm_prev = nlm[tt-1,:]

    dndt_LATROT = sf.dndt_LATROT(nlm_prev, eps,omg) # Here we consider constant large-scale velocity gradients, **but if for any practical time-varying scenario dndt_ij should be calculated inside the loop below!**
    dndt_REG    = sf.dndt_REG(nlm_prev, eps)
    dndt        = dndt_LATROT + dndt_REG

    nlm[tt,:] = nlm_prev + dt*np.matmul(dndt, nlm_prev)

#----------------------
# Aux state vars
#----------------------

# Calculate eigenvalues, principal directions, and enhancement-factors

vecdim = (Nt,3)
eigvals = np.zeros(vecdim, dtype=np.float64)

Eeiej_lin,Eeiej_nlin = np.zeros((Nt,3,3), dtype=np.float64), np.zeros((Nt,3,3), dtype=np.float64)
Epipj_lin,Epipj_nlin = np.zeros((Nt,3,3), dtype=np.float64), np.zeros((Nt,3,3), dtype=np.float64)

e1,e2,e3 = np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64)
p1,p2,p3 = np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64)

for tt in np.arange(0,Nt):

    c = nlm[tt,:]
    
    e1[tt,:],e2[tt,:],e3[tt,:], eigvals[tt,:] = sf.frame(c, 'e')
    p1[tt,:],p2[tt,:],p3[tt,:], _             = sf.frame(c, 'p')

    # Linear (n'=1) mixed Taylor--Sachs enhancements            
    nprime = 1
    Eeiej_lin[tt,:,:] = np.transpose(sf.Eeiej(c, e1[tt,:],e2[tt,:],e3[tt,:], Ecc_lin, Eca_lin, alpha_lin, nprime))
    Epipj_lin[tt,:,:] = np.transpose(sf.Eeiej(c, p1[tt,:],p2[tt,:],p3[tt,:], Ecc_lin, Eca_lin, alpha_lin, nprime))
    
    # Nonlinear (n'=3) Sachs enhancements
    nprime = 3
    Eeiej_nlin[tt,:,:] = np.transpose(sf.Eeiej(c, e1[tt,:],e2[tt,:],e3[tt,:], Ecc_nlin, Eca_nlin, alpha_nlin, nprime))
    Epipj_nlin[tt,:,:] = np.transpose(sf.Eeiej(c, p1[tt,:],p2[tt,:],p3[tt,:], Ecc_nlin, Eca_nlin, alpha_nlin, nprime))
    
#----------------------
# Save
#----------------------

# Save the solution to a netCDF file.
# You can plot the results stored in the netCDF file using "plot.py".

fname = 'solutions/LATROT_%s.nc'%(arg_exp)
ncfile = Dataset(fname,mode='w',format='NETCDF3_CLASSIC') 

# Config
ncfile.tsteps, ncfile.dt, ncfile.L = Nt, dt, L
ncfile.Ecc_lin,  ncfile.Eca_lin,  ncfile.alpha_lin  = Ecc_lin,  Eca_lin,  alpha_lin
ncfile.Ecc_nlin, ncfile.Eca_nlin, ncfile.alpha_nlin = Ecc_nlin, Eca_nlin, alpha_nlin

ncfile.ugrad = ugrad.flatten()

# Dimensions
c_did       = ncfile.createDimension('DOF',     nlm_len)
time_did    = ncfile.createDimension('tstep',   Nt)
eig_did     = ncfile.createDimension('eigval',  3)
dim_did     = ncfile.createDimension('dim',     3)
pair_did    = ncfile.createDimension('pair',    2)

# Variables
myint = np.int32
myflt = np.float64
dimarr_vec = ('lm','lm','lmdyn')

f_lm   = ncfile.createVariable('lm', myint, ('DOF','pair')) 
f_c_re = ncfile.createVariable('c_re', myflt, ('tstep','DOF')) 
f_c_im = ncfile.createVariable('c_im', myflt, ('tstep','DOF')) 
mkvec = lambda field: ncfile.createVariable(field, myflt, ('tstep','dim'))
f_eigvals = ncfile.createVariable('eigvals', myflt, ('tstep','eigval'))

f_Eeiej_lin  = ncfile.createVariable('Eeiej_lin',  myflt, ('tstep','dim','dim'))
f_Epipj_lin  = ncfile.createVariable('Epipj_lin',  myflt, ('tstep','dim','dim'))
f_Eeiej_nlin = ncfile.createVariable('Eeiej_nlin', myflt, ('tstep','dim','dim'))
f_Epipj_nlin = ncfile.createVariable('Epipj_nlin', myflt, ('tstep','dim','dim'))

f_e1, f_e2, f_e3 = mkvec('e1'),  mkvec('e2'),  mkvec('e3')
f_p1, f_p2, f_p3 = mkvec('p1'),  mkvec('p2'),  mkvec('p3')

f_lm[:,:] = lm.T
f_c_re[:,:], f_c_im[:,:] = np.real(nlm), np.imag(nlm)
f_eigvals[:,:] = eigvals
f_Eeiej_lin[:,:,:]  = Eeiej_lin
f_Epipj_lin[:,:,:]  = Epipj_lin
f_Eeiej_nlin[:,:,:] = Eeiej_nlin
f_Epipj_nlin[:,:,:] = Epipj_nlin
f_e1[:,:], f_e2[:,:], f_e3[:,:] = e1, e2, e3
f_p1[:,:], f_p2[:,:], f_p3[:,:] = p1, p2, p3

ncfile.close(); 
print('Solution dumped in %s'%fname)
print('Plot result:\npython3 plot-fabric-evolution-latrot.py %s'%(arg_exp))

