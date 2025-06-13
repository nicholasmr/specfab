import numpy as np
import matplotlib.pyplot as plt
from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
lm, nlm_len = sf.init(4)

### Determine deformation gradient F

# Pure shear 
axis   = 2 # axis of compression/extension (0=x, 1=y, 2=z)
q      = 0 # deformation asymmetry
tau_ps = 1/np.log(2) # e-folding time scale
t_ps   = 1 # time at which deformed parcel is sought
F_ps   = sf.pureshear_F(axis, q, tau_ps, t_ps) # deformation gradient tensor 

# Simple shear
plane   = 1 # shear plane (0=yz, 1=xz, 2=xy)
tau_ss  = 1 # characteristic time taken to reach shear strain 45 deg.
t_ss    = 1 # time at which deformed parcel is sought
F_ss    = sf.simpleshear_F(plane, tau_ss, t_ss) # deformation gradient tensor

### Plot
fig = plt.figure(figsize=(6,6))
ax1 = plt.subplot(121, projection='3d')
ax2 = plt.subplot(122, projection='3d')
sfplt.plotparcel(ax1, F_ps, azim=35, axscale=1.7, axislabels=True, drawinit=True)
sfplt.plotparcel(ax2, F_ss, azim=35, axscale=1.7, axislabels=True, drawinit=True)
ax1.set_title(r'$\epsilon_{zz}=%.2f$'%(sf.F_to_strain(F_ps)[2,2]))
ax2.set_title(r'$\gamma=%.0f$ deg.'%(np.rad2deg(sf.simpleshear_gamma(tau_ss, t_ss))))
plt.savefig('deformed-parcel.png', dpi=175, pad_inches=0.1, bbox_inches='tight')
