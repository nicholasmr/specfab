# N. Rathmann, 2025

import numpy as np
from specfabpy import specfab as sf
sf.init(4)
ex,ey,ez = np.eye(3)

"""
Test the *relative* viscosity matrix model of Elmer/Ice GOLF 

The 6x6 etaij matrix is defined so that
    
    StressVec = 2*eta0 * matmul(etaij, StrainRateVec)
    
where 
    
    StressVec     = [Sxx, Syy, Szz, Sxy, Syz, Sxz]
    StrainRateVec = [Dxx, Dyy, Dzz, Dxy, Dyz, Dxz]
    
Notice that for isotropic ice, Glen's law must be recovered, implying: matmul(etaij, StrainRateVec) = StrainRateVec
"""

### isotropic
a2 = np.eye(3)/3 

### single max
#a2     = np.diag([1,0,0]) # z single max
#a2     = np.diag([0,1,0]) # z single max
#a2     = np.diag([0,0,1]) # z single max

### girdle
#a2     = np.diag([0.5,0.5,0]) 
#a2     = np.diag([0.5,0,0.5]) 
#a2     = np.diag([0,0.5,0.5]) 

etaij = sf.etaij_GOLF(a2)
print(etaij)

etaij_diag = np.diag(etaij)
#print(etaij_diag)
eta_xxxx, eta_yyyy, eta_zzzz, eta_xyxy, eta_yzyz, eta_xzxz = etaij_diag # unpack components

"""
Verify that for an isotropic fabric: matmul(etaij, StrainRateVec) = StrainRateVec
"""

print('\n*** Testing a2 = isotropic...')

a2_iso = np.eye(3)/3
etaij_iso = sf.etaij_GOLF(a2_iso)
tovec = lambda M: np.array([M[0,0],M[1,1],M[2,2], M[0,1],M[1,2],M[0,2]]) # vectorize matrix

epsxz = np.einsum('i,j',ex,ez) + np.einsum('i,j',ez,ex)
epsxz_vec = tovec(epsxz)
print(epsxz_vec)
print(np.matmul(etaij_iso,epsxz_vec))

epszz = np.eye(3)/3 - np.einsum('i,j',ez,ez)
epszz_vec = tovec(epszz)
print(epszz_vec)
print(np.matmul(etaij_iso,epszz_vec))

"""
To rotate the frame into curvilinear coordinates, do

    a2 = np.matmul(R,np.matul(a2,R.T))

where R is the rotation matrix that, at each point in space, rotates the cartesian coordinate system into a frame where the "x" direction becomes tangential to flow lines ("t" direction), and the "y" direction becomes normal to flow lines ("n" direction).
"""

