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
# Construct gaunt coefs
#---------------------------------

dl = 2
def lm_generator(L): return [(int(l),int(m-l),int(m)) for l in np.arange(0,L+1,dl) for m in np.arange(0,2*l+1,1)] # NOTE: must be naitive int() vals for gaunt() to behave correctly.
lm     = lm_generator(L) # full l,m list 
#lm_dyn = lm_generator(2) # dynamic catalyst l,m list :: if lattice rotation only
lm_dyn = lm_generator(4) # dynamic catalyst l,m list  :: required for DDRX
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
# Save ".f90 body"
#---------------------------------

def savecoef(f, f90var, gaunt):
    for ii, (li,mi,mii) in enumerate(lm): # i'th Psi coef to be determined (conjugate term in overlap integral)
        for jj, (lj,mj,mjj) in enumerate(lm): # j'th coef in sum over Psi
            for kk, (lk,mk,mkk) in enumerate(lm_dyn): # k'th coef for "catalyst" terms in tripple products
                if np.abs(gaunt[ii,jj,kk])>1e-20:
                    f.write("%s(%i,%i,%i) = %f\n" % (f90var,kk+1,jj+1,ii+1,np.real(gaunt[ii,jj,kk])) )

f = open("gaunt__head.f90", "w")
f.write('! L = %i \n'%(L))
f.write('real(kind=dp), dimension(%i,%i,%i) :: GC=0.0, GCm=0.0, GC_m1=0.0, GC_p1=0.0 \n'%(len(lm_dyn),Nc,Nc))
f.close()

f = open("gaunt__body.f90", "w")
savecoef(f, 'GC', GC)
savecoef(f, 'GCm', GCm)
savecoef(f, 'GC_m1', GC_m1)
savecoef(f, 'GC_p1', GC_p1)
f.close()

