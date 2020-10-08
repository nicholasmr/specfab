#!/usr/bin/python
# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

import sys
import numpy as np
from numpy import sqrt, power, log, sum, tensordot, multiply, conj, rad2deg, deg2rad
from sympy.physics.wigner import gaunt, wigner_3j
from sympy.parsing.sympy_parser import parse_expr
from netCDF4 import Dataset 

L = int(sys.argv[1])

#---------------------------------
# Construct gaunt coef matrix
#---------------------------------

dl = 2
def lm_generator(L): return [(int(l),int(m-l),int(m)) for l in np.arange(0,L+1,dl) for m in np.arange(0,2*l+1,1)] # NOTE: must be naitive int() vals for gaunt() to behave correctly.
lm     = lm_generator(L) # full l,m list 
lm_dyn = lm_generator(2) # dynamic catalyst l,m list 
Nc = len(lm)
            
def gauntmatrix(lm_j, lm_cat, wk=None, wj=None, delta_mj=0):

    GC = np.zeros((Nc, len(lm_j), len(lm_cat)), dtype=np.complex64)
    if wj is None: wj = np.ones(len(lm_j))
    if wk is None: wk = np.ones(len(lm_cat))

    if 1: print("GC matrix with L=%i of shape "%(L),GC.shape, '...')
    lmkk = np.array([ [li,mi] for li,mi,mii in lm_cat ])
    
    for ii, (li,mi,mii) in enumerate(lm): # i'th Psi coef to be determined (conjugate term in overlap integral)
        print('gauntmatrix: outer index %i of %i (l=%i)'%(ii,Nc,li))
        for jj, (lj,mj,mjj) in enumerate(lm_j): # j'th coef in sum over Psi
            for kk, (lk,mk,mkk) in enumerate(lm_cat): # k'th coef for "catalyst" terms in tripple products
                GC[ii,jj,kk] = wj[jj]*wk[kk]*overlap(lk,lj,li, mk,mj+delta_mj,mi)
    
    return GC

def s(l,m): return 0.5*sqrt((l-m)*(l+m+1));

def overlap(lk,lj,li, mk,mj,mi, GAUNT_PRECISION=6): 
    return np.complex64( (-1)**abs(mi) * gaunt(int(lk),int(lj),int(li), int(+mk),int(+mj),int(-mi), prec=GAUNT_PRECISION) ) # Last {l,m} index pair are associated with the conjugate term.

GC      = gauntmatrix(lm, lm_dyn) # for dynamics (rotation)
GCm     = gauntmatrix(lm, lm_dyn, wj=[mj for (lj,mj,mjj) in lm]) # ... m-weighted
GC_p1   = gauntmatrix(lm, lm_dyn, wj=[s(lj,+mj) for (lj,mj,mjj) in lm], delta_mj=+1)
GC_m1   = gauntmatrix(lm, lm_dyn, wj=[s(lj,-mj) for (lj,mj,mjj) in lm], delta_mj=-1)

#---------------------------------
# Save to netCDF
#---------------------------------

fname = 'gauntcoefs_L%02i.nc'%(L)
print('Saving to %s'%(fname))
ncfile = Dataset(fname,mode='w',format='NETCDF3_CLASSIC') 

dim_lmdyn = ncfile.createDimension('lmdyn', len(lm_dyn))    
dim_lm    = ncfile.createDimension('lm',  len(lm))

dtype = np.float32
dimarr = ('lm','lm','lmdyn')
f_GC     = ncfile.createVariable('GC',    dtype, dimarr) 
f_GCm    = ncfile.createVariable('GCm',   dtype, dimarr)
f_GC_p1  = ncfile.createVariable('GC_p1', dtype, dimarr) 
f_GC_m1  = ncfile.createVariable('GC_m1', dtype, dimarr)  

f_GC[:,:,:]    = np.real(GC)
f_GCm[:,:,:]   = np.real(GCm)
f_GC_p1[:,:,:] = np.real(GC_p1)
f_GC_m1[:,:,:] = np.real(GC_m1)

ncfile.close(); 
print('Dumped structs to %s'%fname)

