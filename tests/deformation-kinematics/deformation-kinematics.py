# N. M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Test parcel deformation kinematic routines, including plotting.
Use for specfab documentation site.
"""

import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) 

axis = 2 # axis of shortening (tauc>0) or lengthening (tauc<0): 0=x, 1=y, 2=z
tauc = 1 # time taken in seconds for parcel to reduce to half (50%) height if tauc>0, or abs(time) taken for parcel to double in height (200%) if tauc<0.
q    = 0 # asymmetry parameter for shortening (if tauc>0) or lengthening (if tauc<0)

tau = tauc/np.log(2)                     # corresponding e-folding time
ugrad = sf.pureshear_ugrad(axis, q, tau) # velocity gradient
D, W = sf.ugrad_to_D_and_W(ugrad)        # strain-rate and spin tensor

t = 1 # some specific time of interest
r   = sf.pureshear_r(tau, t)          # scaling parameter "r" at time t
F   = sf.pureshear_F(axis, q, tau, t) # deformation gradient tensor at time t
eps = sf.F_to_strain(F)               # strain tensor at time t

print(r)
print(eps)

######################

import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) 

plane = 1 # plane of shear: 0=yz, 1=xz, 2=xy
tau   = 1 # time taken in seconds for parcel to a reach shear strain of 1 (45 deg shear angle)

ugrad = sf.simpleshear_ugrad(plane, tau) # velocity gradient
D, W = sf.ugrad_to_D_and_W(ugrad)        # strain-rate and spin tensor

t = 1 # some specific time of interest
gamma = sf.simpleshear_gamma(tau, t)    # shear angle at time t
F     = sf.simpleshear_F(plane, tau, t) # deformation gradient tensor at time t
eps   = sf.F_to_strain(F)               # strain tensor at time t

print(np.rad2deg(gamma))
print(eps)

########################

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

#########################


import numpy as np
import matplotlib.pyplot as plt
from specfabpy import specfab as sf
from specfabpy import plotting as sfplt

lm, nlm_len = sf.init(6)

### CPO to be plotted

a2 = np.diag([0,0,1]) # CPO characterized by a^(2)
nlm = sf.a2_to_nlm(a2) # vector of expansion coefficients

### Setup axes and projection

geo, prj = sfplt.getprojection(rotation=45, inclination=45)
fig = plt.figure(figsize=(2,2))
ax = plt.subplot(111, projection=prj)
ax.set_global() # ensure entire S^2 is shown

### Plot

lvlset = 'iso-up' # default level set: lowest tick/level is the value of an isotropic distribution
lvlset = (np.linspace(0,0.8,9), lambda x,p:'%.1f'%x) # custom level set: (list of levels, how to format colorbar tick labels)
sfplt.plotODF(nlm, lm, ax, cmap='Greys', cblabel='ODF', lvlset=lvlset) # plot distribution (see src/specfabpy/plotting.py for API)
sfplt.plotcoordaxes(ax, geo, color=sfplt.c_dred) # plot coordinate axes (see src/specfabpy/plotting.py for API)

plt.savefig('ODF-plot.png', dpi=175, pad_inches=0.1, bbox_inches='tight')

