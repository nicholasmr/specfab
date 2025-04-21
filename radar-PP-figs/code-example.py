import numpy as np
from scipy.spatial.transform import Rotation
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) # L=4 is sufficient here

### Determine <c_i c_j> from radar-derived Delta lambda
l1 = 0   # lambda_1 = 0 (Gerber's approximation)
dl = 0.5 # Delta lambda = lambda_2 - lambda_1
a2 = np.diag([l1, l1+dl, 1-dl-2*l1]) # second-order structure tensor, <c_i c_j>, in eigenframe
m1, m2, z = np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1]) # eigenvectors

### Rotate <c_i c_j> into a rotationally-symmetric frame about z
Rm1 = Rotation.from_rotvec(np.pi/2 * m1).as_matrix() # Rotate 90 deg about m1 eigenvector
Rm2 = Rotation.from_rotvec(np.pi/2 * m2).as_matrix() # Rotate 90 deg.about m2 eigenvector
if dl < 0.4:         a2_vs = a2 # Already in rotationally-symmetric frame about z
if 0.4 <= dl <= 0.6: a2_vs = np.matmul(Rm2,np.matmul(a2,Rm2.T)) # Rotate vertical (m2--z) girdle into horizontal (m1--m2) girdle
if dl > 0.6:         a2_vs = np.matmul(Rm1,np.matmul(a2,Rm1.T)) # Rotate horizontal (m2) single-maximum into vertical (z) single-maximum

### Determine \hat{n}_4^0 (= n_4^0/n_0^0) from \hat{n}_2^0 (= n_2^0/n_0^0) in rotationally-symmetric frame about z
nhat20 = (a2_vs[2,2]- 1/3)/(2/15*np.sqrt(5)) # azz -> nhat20
nhat40 = sf.nhat40_empcorr_ice(nhat20)[0]

### Construct nlm (spectral CPO state vector) in rotationally-symmetric frame about z
nlm_vs = np.zeros(nlm_len, dtype=np.complex128) 
n00 = 1/np.sqrt(4*np.pi) # only grain-number normalized distribution is known, so must integrate to 1 over S^2.
nlm_vs[0]  = n00
nlm_vs[3]  = nhat20*n00
nlm_vs[10] = nhat40*n00

### Rotate spectral CPO state back to origional (m1,m2,z) eigenframe 
if dl < 0.4:         nlm = nlm_vs # Already in vertical symmetric frame
if 0.4 <= dl <= 0.6: nlm = sf.rotate_nlm(nlm_vs, -np.pi/2, 0) # Rotate horizontal (m1--m2) girdle back into vertical (m2--z) girdle
if dl > 0.6:         nlm = sf.rotate_nlm(sf.rotate_nlm(nlm_vs, -np.pi/2, 0), 0 ,-np.pi/2) # Rotate vertical (z) single-maximum back into horizontal (m2) single-maximum

### Calculate eigenenhancements
# Transversely isotropic monocrystal parameters for ice (Rathmann & Lilien, 2021)
n_grain   = 1        # Power-law exponent: n=1 => linear grain rheology, nonlinear (n>1) is unsupported
Eij_grain = (1, 1e3) # Grain eigenenhancements (Ecc,Eca) for compression along c-axis (Ecc) and for shear parallel to basal plane (Eca)
alpha     = 0.0125   # Taylor--Sachs weight
# Tuple of eigenenhancements (bulk enhancement factors w.r.t. m1, m2, z)
e1, e2, e3 = m1, m2, z
Eij = sf.Eij_tranisotropic(nlm, e1,e2,e3, Eij_grain,alpha,n_grain) # (E_{m1,m1},E_{m2,m2},E_{zz},E_{m2,z),E_{m1,z},E_{m1,m2})
# To calculate bulk enhancement factors w.r.t. other axes of deformation/stress, change (e1,e2,e3) accordingly.


#-------------------------------
#-------------------------------

import copy, sys, code # code.interact(local=locals())
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

geo, prj = sfplt.getprojection(rotation=55, inclination=45)

scale = 2.5
fig = plt.figure(figsize=(3.4/2*1.45*scale,0.86*scale))
gs = fig.add_gridspec(1,4)
al = 0.04
ar = 0.02
gs.update(left=al, right=1-ar, top=1.00, bottom=0.22, wspace=0.4, hspace=0.4)
ax1 = fig.add_subplot(gs[0,0], projection=prj); ax1.set_global(); 
ax2 = fig.add_subplot(gs[0,1], projection=prj); ax2.set_global(); 
ax3 = fig.add_subplot(gs[0,2], projection=prj); ax3.set_global(); 
ax4 = fig.add_subplot(gs[0,3], projection=prj); ax4.set_global(); 

isgirdle = 0.4 <= dl <= 0.6
#lvlset = 'iso-up'
lvlset = [np.linspace(0.15,0.25,7) if isgirdle else np.linspace(0.15,0.65,7), lambda x,p:'%.2f'%x]
cbtickintvl = 3 if isgirdle else 3
kw = dict(lvlset=lvlset, cbtickintvl=cbtickintvl, cblabel='ODF')
sfplt.plotODF(sf.a2_to_nlm(a2),lm,    ax1, **kw)
sfplt.plotODF(sf.a2_to_nlm(a2_vs),lm, ax2, **kw)
sfplt.plotODF(nlm_vs,lm,              ax3, **kw)
sfplt.plotODF(nlm,lm,                 ax4, **kw)
for axi in [ax1,ax2,ax3,ax4]: sfplt.plotcoordaxes(axi, geo, axislabels=[r'$\vb{m}_1$',r'$\vb{m}_2$',r'$\vb{z}$'])

ax1.set_title(r'\texttt{sf.a2\_to\_nlm(a2)}', fontsize=FS)
ax2.set_title(r'\texttt{sf.a2\_to\_nlm(a2\_vs)}', fontsize=FS)
ax3.set_title(r'\texttt{nlm\_vs}', fontsize=FS)
ax4.set_title(r'\texttt{nlm}', fontsize=FS)

print('Eij = ',Eij)
plt.savefig('code-example-output-dl%.1f.png'%(dl), transparent=True, dpi=140)

