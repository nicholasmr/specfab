# N. Rathmann, 2025

"""
Test Elmer/Ice GOLF viscosity model

The 6x6 viscosity matrix etaij, dimensionless so that
    
    StressVec = 1/B * matmul(etaij, StrainRateVec)
    
where stress and strainrate vectors are the ordered list of components
    
    [xx, yy, zz, xy, yz, xz]
"""

import numpy as np
from specfabpy import specfab as sf
sf.init(4)
ex,ey,ez = np.eye(3)

### isotropic
a2 = np.eye(3)/3 

### single max
a2     = np.diag([1,0,0]) # z single max
#a2     = np.diag([0,1,0]) # z single max
#a2     = np.diag([0,0,1]) # z single max

### girdle
a2     = np.diag([0.5,0.5,0]) 
#a2     = np.diag([0.5,0,0.5]) 
#a2     = np.diag([0,0.5,0.5]) 

etaij = sf.etaij_GOLF(a2)
print(etaij)

etaij_diag = np.diag(etaij)
print(etaij_diag)

eta_xxxx, eta_yyyy, eta_zzzz, eta_xyxy, eta_yzyz, eta_xzxz = etaij_diag # unpack components

"""
To rotate the frame into curvilinear coordinates, do

    a2 = np.matmul(R,np.matul(a2,R.T))

where R is the rotation matrix that, at each point in space, rotates the cartesian coordinate system into a frame where the "x" direction becomes tangential to flow lines ("t" direction), and the "y" direction becomes normal to flow lines ("n" direction).
"""

