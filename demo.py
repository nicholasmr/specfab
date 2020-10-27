# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp
from netCDF4 import Dataset
from specfabpy import specfabpy as sf 

#### NOTICE ####
# To use "specfabpy", you must first compile the Python module by running "make demopy" (compiles the F2PY wrapper for "specfab.f90")
################

#----------------------
# Input arguments
#----------------------

if len(sys.argv) != 4:
    print('usage: python3 %s L nprime ugrad'%(sys.argv[0]))
    sys.exit(0)
    
# Arguments
arg_L      = int(sys.argv[1])   # Series truncation
arg_nprime = int(sys.argv[2])   # Grain rheology exponent
arg_exp    =     sys.argv[3]    # Velocity gradient experiment

L = arg_L
if not(2 <= L <= 40):
    print('L must be 2 <= L <= 40')
    sys.exit(0)

#----------------------
# Grain rheology
#----------------------

# Relevant only for calculating the enhancement-factors resulting from the fabric.

if arg_nprime == 1: nprime, Ecc, Eca = 1, 1, 1e4
if arg_nprime == 3: nprime, Ecc, Eca = 3, 1, 1e2

#----------------------
# Numerics
#----------------------

Nt = 300  # Number of integration points
dt = 0.02 # Time-step size

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

nlm_len = sf.init(L) # nlm_len is the number of fabric expansion coefficients (degrees of freedom).
nlm = np.zeros((Nt,nlm_len), dtype=np.complex128) # The expansion coefficients
lm = sf.get_lm(nlm_len) # The (l,m) values corresponding to the coefficients in "nlm".
nu = sf.get_nu(eps) # Diffusion coefficient used for regularization.

# Initial fabric state
if arg_exp[0:2] == 'rr':
    # Init with some anisotropy
    nlm[0,0] = np.sqrt(1/2 + 0j)
    nlm[0,2] = np.sqrt(1/2 + 0j)
else:
    # Init with isotropy
    nlm[0,0] = 1/np.sqrt(4*np.pi) # Normalized such that N(t=0) = 1

#----------------------
# Integrate
#----------------------

print('Numerics: Nt=%i, dt=%f, L=%i (nlm_len=%i)'%(Nt,dt,L,nlm_len))

# Euler scheme
for tt in np.arange(1,Nt):
    nlm[tt,:] = nlm[tt-1,:] + dt*np.tensordot(sf.get_dndt_ij(nlm_len, eps, omg), nlm[tt-1,:], axes=(-1,0))

#----------------------
# Aux state vars
#----------------------

# Calculate eigenvalues, principal directions, and enhancement-factors

vecdim = (Nt,3)
eigvals = np.zeros(vecdim, dtype=np.float64)
Eeiej   = np.zeros((Nt,3,3), dtype=np.float64)
Epijqij = np.zeros(vecdim, dtype=np.float64)
e1,e2,e3    = np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64)
p23,p12,p13 = np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64)
q23,q12,q13 = np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64)

for tt in np.arange(0,Nt):
    c = nlm[tt,:]
    
    e1[tt,:],e2[tt,:],e3[tt,:], eigvals[tt,:] = sf.get_eigenframe(c)
    Eeiej[tt,:,:] = np.transpose(sf.get_eeiej(c, Ecc,Eca,nprime))
    
    p23[tt,:],p12[tt,:],p13[tt,:], q23[tt,:],q12[tt,:],q13[tt,:] = sf.get_pqframe(c)
    Epijqij[tt,:] = sf.get_epijqij(c, Ecc,Eca,nprime)
    
#----------------------
# Save
#----------------------

# Save the solution to a netCDF file.
# You can plot the results stored in the netCDF file using "plot.py".

fname = 'solutions/solution_n%i_%s.nc'%(nprime, arg_exp)
ncfile = Dataset(fname,mode='w',format='NETCDF3_CLASSIC') 

# Config
ncfile.tsteps, ncfile.dt, ncfile.nu, ncfile.L = Nt, dt, nu, L
ncfile.nprime, ncfile.Ecc, ncfile.Eca = nprime, Ecc, Eca
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
f_Eeiej   = ncfile.createVariable('Eeiej',   myflt, ('tstep','dim','dim'))
f_Epijqij = ncfile.createVariable('Epijqij', myflt, ('tstep','dim'))
f_e1, f_e2, f_e3    = mkvec('e1'),  mkvec('e2'),  mkvec('e3')
f_p23, f_p12, f_p13 = mkvec('p23'), mkvec('p12'), mkvec('p13')
f_q23, f_q12, f_q13 = mkvec('q23'), mkvec('q12'), mkvec('q13')

f_lm[:,:] = lm.T
f_c_re[:,:], f_c_im[:,:] = np.real(nlm), np.imag(nlm)
f_eigvals[:,:] = eigvals
f_Eeiej[:,:,:] = Eeiej
f_Epijqij[:,:] = Epijqij
f_e1[:,:], f_e2[:,:], f_e3[:,:]    = e1, e2, e3
f_p23[:,:], f_p12[:,:], f_p13[:,:] = p23, p12, p13
f_q23[:,:], f_q12[:,:], f_q13[:,:] = q23, q12, q13

ncfile.close(); 
print('Solution dumped in %s'%fname)
print('Plot result:\npython3 plot.py %i %s'%(nprime, arg_exp))

