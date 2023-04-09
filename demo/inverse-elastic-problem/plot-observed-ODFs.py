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
from plottools import *

### Init figure

inclination = 50 # view angle
rot = -40 # view angle
#rot, inclination = -90*1, 0 # debug, view trom top--down
prj = ccrs.Orthographic(rot, 90-inclination)
geo = ccrs.Geodetic()

scale = 0.45
fig = plt.figure(figsize=(10*scale,5*scale))
gs = gridspec.GridSpec(1, 3) 
a=0.025
gs.update(top=0.99, bottom=0.22, left=a, right=1-a, wspace=0.35)

ax = [fig.add_subplot(gs[0, jj], projection=prj) for jj in range(3)]

### Plot

Lutz = Lutz_etal_2022(verbose=False)

for ii, exprnum in enumerate([3,7,10]):

    (nlm, qlat,qlon) = Lutz.get_nlm(exprnum) 

    ax[ii].set_global() 

    tickintvl, lvls = 3, np.linspace(0.0,0.6,7)
    kwargs_ODF = {'cmap':'Greys', 'cblabel':'$\psi/N$', 'cbaspect':8.5, 'cbfrac':0.065, 'lvls':lvls, 'tickintvl':tickintvl}
    plot_ODF(nlm, lm_L4, ax=ax[ii], **kwargs_ODF)

    qlatd, qlond = get_deg(qlat, qlon)
    ax[ii].plot(qlond, qlatd, ls='none', marker='o', markersize=0.45, c='#33a02c', transform=geo) 
    plot_unitaxes(ax[ii], geo)
    ax[ii].set_title(r'sample %03i'%(exprnum), pad=10, fontsize=FS)

fname = 'plots/observed-ODFs.png'
#fname = 'plots/observed-ODFs.pdf'
print('** Saving %s'%(fname))
plt.savefig(fname, dpi=250)  
  
