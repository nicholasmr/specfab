# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

""" 
Plot observed ODFs of Thomas et al. (2021) samples #003, #007, and #010
"""

import copy, sys, code # code.interact(local=locals())
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from experiments import *
from inverseproblem import lm_L4

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

### Init figure

scale = 0.45
fig = plt.figure(figsize=(10*scale,5*scale))
gs = gridspec.GridSpec(1, 3) 
a=0.025
gs.update(top=0.99, bottom=0.22, left=a, right=1-a, wspace=0.35)

geo, prj = sfplt.getprojection(rotation=50-90, inclination=50)
ax = [fig.add_subplot(gs[0, jj], projection=prj) for jj in range(3)]

### Plot

Lutz = Lutz_etal_2022(verbose=False)

for ii, exprnum in enumerate([3,7,10]):

    (nlm, qlat,qlon) = Lutz.get_nlm(exprnum) 

    ax[ii].set_global() # show entire S^2

    lvlset = (np.linspace(0.0,0.6,7), lambda x,p:'%.1f'%x)
    sfplt.plotODF(nlm, lm_L4, ax[ii], lvlset=lvlset, cbtickintvl=3, cblabel='$\psi/N$', cbaspect=8.5, cbfraction=0.065)

    ax[ii].plot(np.rad2deg(qlon), np.rad2deg(qlat), ls='none', marker='o', markersize=0.45, c='#33a02c', transform=geo) 

    sfplt.plotcoordaxes(ax[ii], geo, axislabels='vuxi', color=sfplt.c_dred)  

    ax[ii].set_title(r'sample %03i'%(exprnum), pad=10, fontsize=FS)

fname = 'plots/observed-ODFs.png'
print('** Saving %s'%(fname))
plt.savefig(fname, dpi=250)  
  
