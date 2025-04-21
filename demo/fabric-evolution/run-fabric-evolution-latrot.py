# N. M. Rathmann <rathmann@nbi.ku.dk>, 2020-2023

import sys, os, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp
from netCDF4 import Dataset

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import constants as sfconst

#----------------------
# Experiment selection
#----------------------

if len(sys.argv) != 2:
    print('usage: python3 %s ugrad'%(sys.argv[0]))
    sys.exit(0)
    
exp     = sys.argv[1]
exptype = exp[0:2]
ijstr   = exp[-2:]

ii_to_axis  = lambda ij: 0 if (ij=='xx') else (1 if (ij=='yy') else 2)
ij_to_plane = lambda ij: 0 if (ij=='yz') else (1 if (ij=='xz') else 2)

if exptype == 'ue': # Uniaxial extension
    DK = dict(type='ps', tau=-1, q=0,  axis=ii_to_axis(ijstr))
    strain_target = 3

if exptype == 'uc': # Uniaxial compression
    DK = dict(type='ps', tau=1, q=0,  axis=ii_to_axis(ijstr))
    strain_target = -0.98

if exptype == 'cc': # Confined compression
    DK = dict(type='ps', tau=1, q=+1, axis=ii_to_axis(ijstr))
    strain_target = -0.98

if exptype == 'ss': # Simple shear
    DK = dict(type='ss', tau=1, plane=ij_to_plane(ijstr))
    strain_target = np.deg2rad(70)
    
if exptype == 'rr': # Ridgid rotation
    DK = dict(type='rr', tau=1, plane=ij_to_plane(ijstr))
    strain_target = np.deg2rad(90)

#----------------------
# Model integration
#----------------------

Nt = 50 # Number of time steps
L  = 8 # Spectral truncation

lm, nlm_len = sf.init(L)
nlm, F, time, ugrad = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, iota=1, nu=1)

#----------------------
# Determine eigenvalues, principal directions, and enhancement factors
#----------------------

# Empty structure to fill
vecdim = (Nt+1,3)
eigvals  = np.zeros(vecdim)
m1,m2,m3 = np.zeros(vecdim),np.zeros(vecdim),np.zeros(vecdim)
p1,p2,p3 = np.zeros(vecdim),np.zeros(vecdim),np.zeros(vecdim)
vecdim = (Nt+1,6)
Eij_lin, Eij_nlin   = np.zeros(vecdim), np.zeros(vecdim)
Epij_lin, Epij_nlin = np.zeros(vecdim), np.zeros(vecdim)

# Grain parameters
(Eij_grain_lin,  alpha_lin,  n_grain_lin)  = sfconst.ice['viscoplastic']['linear']
(Eij_grain_nlin, alpha_nlin, n_grain_nlin) = sfconst.ice['viscoplastic']['nonlinear']

viscparams = (Eij_grain_lin, alpha_lin, n_grain_lin)

for tt in np.arange(0,Nt+1):

    c = nlm[tt,:]
    
    m1[tt,:],m2[tt,:],m3[tt,:], eigvals[tt,:] = sf.frame(c, 'e')
    p1[tt,:],p2[tt,:],p3[tt,:], _             = sf.frame(c, 'p')

    # Linear (n'=1) mixed Taylor--Sachs enhancements            
    Eij_lin[tt,:]  = sf.Eij_tranisotropic(c, m1[tt,:],m2[tt,:],m3[tt,:], Eij_grain_lin, alpha_lin, n_grain_lin)
    Epij_lin[tt,:] = sf.Eij_tranisotropic(c, p1[tt,:],p2[tt,:],p3[tt,:], Eij_grain_lin, alpha_lin, n_grain_lin)
    
    # Nonlinear (n'=3) Sachs enhancements
    Eij_nlin[tt,:]  = sf.Eij_tranisotropic(c, m1[tt,:],m2[tt,:],m3[tt,:], Eij_grain_nlin, alpha_nlin, n_grain_nlin)
    Epij_nlin[tt,:] = sf.Eij_tranisotropic(c, p1[tt,:],p2[tt,:],p3[tt,:], Eij_grain_nlin, alpha_nlin, n_grain_nlin)
    
#----------------------
# Save
#----------------------

fname = 'solutions/LROT_%s.nc'%(exp)
ncfile = Dataset(fname,mode='w',format='NETCDF3_CLASSIC') 

# Config
dt = time[1]-time[0]
ncfile.tsteps, ncfile.dt, ncfile.L = Nt, dt, L
ncfile.Ecc_lin,  ncfile.Eca_lin,  ncfile.alpha_lin  = Eij_grain_lin[0],  Eij_grain_lin[1],  alpha_lin
ncfile.Ecc_nlin, ncfile.Eca_nlin, ncfile.alpha_nlin = Eij_grain_nlin[0], Eij_grain_nlin[1], alpha_nlin

ncfile.ugrad = ugrad.flatten()

# Dimensions
c_did       = ncfile.createDimension('DOF',     nlm_len)
time_did    = ncfile.createDimension('tstep',   Nt+1)
eig_did     = ncfile.createDimension('eigval',  3)
dim_did     = ncfile.createDimension('dim',     3)
dim6_did    = ncfile.createDimension('dim6',    6)
pair_did    = ncfile.createDimension('pair',    2)

# Variables
myint = np.int32
myflt = np.float64
dimarr_vec = ('lm','lm','lmdyn')

f_lm   = ncfile.createVariable('lm', myint, ('DOF','pair')) 
f_c_re = ncfile.createVariable('c_re', myflt, ('tstep','DOF')) 
f_c_im = ncfile.createVariable('c_im', myflt, ('tstep','DOF')) 
f_eigvals = ncfile.createVariable('eigvals', myflt, ('tstep','eigval'))

f_Eij_lin   = ncfile.createVariable('Eij_lin',   myflt, ('tstep','dim6'))
f_Epij_lin  = ncfile.createVariable('Epij_lin',  myflt, ('tstep','dim6'))
f_Eij_nlin  = ncfile.createVariable('Eij_nlin',  myflt, ('tstep','dim6'))
f_Epij_nlin = ncfile.createVariable('Epij_nlin', myflt, ('tstep','dim6'))

mkvec = lambda field: ncfile.createVariable(field, myflt, ('tstep','dim'))
f_m1, f_m2, f_m3 = mkvec('m1'),  mkvec('m2'),  mkvec('m3')
f_p1, f_p2, f_p3 = mkvec('p1'),  mkvec('p2'),  mkvec('p3')

f_lm[:,:], f_c_re[:,:], f_c_im[:,:] = lm.T, np.real(nlm), np.imag(nlm)
f_Eij_lin[:,:],  f_Eij_nlin[:,:]  = Eij_lin,  Eij_nlin
f_Epij_lin[:,:], f_Epij_nlin[:,:] = Epij_lin, Epij_nlin
f_m1[:,:], f_m2[:,:], f_m3[:,:] = m1, m2, m3
f_p1[:,:], f_p2[:,:], f_p3[:,:] = p1, p2, p3
f_eigvals[:,:] = eigvals

ncfile.close(); 
print('Solution dumped in %s'%fname)
print('Plot result: python3 plot-fabric-evolution-latrot.py %s'%(exp))

