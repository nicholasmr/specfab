# N. Rathmann, 2025

import numpy as np
from specfabpy import specfab as sf
sf.init(4)

"""
Test the *relative* viscosity matrix model of Elmer/Ice GOLF 

The 6x6 matrix "etaij" is defined so that
    
    Svec = 2*eta0 * matmul(etaij, Dvec)
    
where Svec and Dvec are the vectorized deviatoric stress and strain rate tensors:
    
    Svec = [Sxx, Syy, Szz, Sxy, Syz, Sxz]
    Dvec = [Dxx, Dyy, Dzz, Dxy, Dyz, Dxz]
    
Notice that for isotropic ice, Glen's law must be recovered, implying:  matmul(etaij, Dvec) = Dvec
"""

ex,ey,ez = np.eye(3) # Cartesian basis vectors
vectorize = lambda M: np.array([M[0,0],M[1,1],M[2,2], M[0,1],M[1,2],M[0,2]]) # vectorize matrix

# Ideal fabric states

a2_iso  = np.eye(3)/3 
a2_smax = lambda mi: np.einsum('i,j',mi,mi) # single max with rot. symmetry axis mi
a2_gdl  = lambda mi: (np.eye(3)-a2_smax(mi))/2 # girdle with rot. symmetry axis mi

# Print etaij for some test case

a2_test = a2_smax(ez)

etaij = sf.etaij_GOLF(a2_test)
etaij_diag = np.diag(etaij)
eta_xxxx, eta_yyyy, eta_zzzz, eta_xyxy, eta_yzyz, eta_xzxz = etaij_diag # unpack components
print(etaij)
#print(etaij_diag)

"""
Verify that for an isotropic fabric: matmul(etaij, Dvec) = Dvec
"""

print('\n*** Testing a2 = isotropic...')

a2_iso = np.eye(3)/3
etaij_iso = sf.etaij_GOLF(a2_iso)

epsxz = np.einsum('i,j',ex,ez) + np.einsum('i,j',ez,ex)
epsxz_vec = vectorize(epsxz)
print(epsxz_vec)
print(np.matmul(etaij_iso,epsxz_vec))

epszz = np.eye(3)/3 - np.einsum('i,j',ez,ez)
epszz_vec = vectorize(epszz)
print(epszz_vec)
print(np.matmul(etaij_iso,epszz_vec))

"""
To rotate the frame into curvilinear coordinates, do

    a2 = np.matmul(R,np.matul(a2,R.T))

where R is the rotation matrix that, at each point in space, rotates the cartesian coordinate system into a frame where the "x" direction becomes tangential to flow lines ("t" direction), and the "y" direction becomes normal to flow lines ("n" direction).
"""

