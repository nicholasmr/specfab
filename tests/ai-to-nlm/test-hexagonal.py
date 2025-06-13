#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2025-

"""
Verify that constructing <a^6> from discrete ensemble of a-axes 
is sufficient for capturing hexagonal symmetry if nlm = sf.a6_to_nlm() is used.
"""

import copy, sys, os, code # code.interact(local=locals())
import numpy as np
from scipy.spatial.transform import Rotation as R

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
import cartopy.crs as ccrs

FS = sfplt.setfont_tex(fontsize=12)
lm, nlm_len = sf.init(8)

ex,ey,ez = np.eye(3)

a2 = np.zeros((3,)*2)
a4 = np.zeros((3,)*4)
a6 = np.zeros((3,)*6)

I = np.arange(3)
N = len(I)
for ii in I:
    ang = ii*np.deg2rad(60)
    print('a-axis angle: %.1f'%(np.rad2deg(ang)))
    r = R.from_rotvec(ang*ez)
    a_ = np.matmul(r.as_matrix(), ex)
    a2 += np.einsum('i,j',         *[a_,]*2)/N
    a4 += np.einsum('i,j,k,l',     *[a_,]*4)/N
    a6 += np.einsum('i,j,k,l,m,n', *[a_,]*6)/N

nlm_2 = np.zeros(nlm_len, dtype=np.complex128)
nlm_4 = np.zeros(nlm_len, dtype=np.complex128)
nlm_6 = np.zeros(nlm_len, dtype=np.complex128)

nlm_2[:sf.L2len] = sf.a2_to_nlm(a2) 
nlm_4[:sf.L4len] = sf.a4_to_nlm(a4) 
nlm_6[:sf.L6len] = sf.a6_to_nlm(a6) 

### Plot

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.gridspec as gridspec

geo, prj = sfplt.getprojection(rotation=45-10, inclination=50)
prj2 = ccrs.AlbersEqualArea(central_longitude=0) # use another projection for a axes

fig = plt.figure(figsize=(5,3))
gs = gridspec.GridSpec(2, 3)
gs.update(left=0.07, right=0.97, top=0.96, bottom=0.1, hspace=0.4)
axi = [fig.add_subplot(gs[0,ii], projection=prj ) for ii in range(3)]
axj = [fig.add_subplot(gs[1,ii], projection=prj2) for ii in range(3)]

for axset in (axi,axj):
    for ax in axset: ax.set_global()

lvlset = (np.linspace(0,0.8,9), lambda x,p:'%.1f'%x) 

for ii,nlm in enumerate((nlm_2,nlm_4,nlm_6)):
           
    sfplt.plotODF(nlm, lm, axi[ii], cblabel=r'ODF', lvlset=lvlset)
    sfplt.plotcoordaxes(axi[ii], geo, axislabels='vuxi', color='tab:red')

    sfplt.plotODF(nlm, lm, axj[ii], cblabel=r'ODF', lvlset=lvlset)
    sfplt.plotcoordaxes(axj[ii], geo, axislabels='vuxi', color='tab:red')
   
    axi[ii].set_title(r'\texttt{a%i_to_nlm()}'%((ii+1)*2))
   
plt.savefig('test-hexagonal.png', dpi=150, bbox_inches='tight', pad_inches=0.05)

