# Nicholas Rathmann, 2025

import numpy as np
import code # code.interact(local=locals())

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage[sc]{mathpazo} \usepackage{siunitx} \DeclareSIUnit\year{a}'

lm, nlm_len = sf.init(4)

scale = 1.5
fig = plt.figure(figsize=(3.4*scale,1.1*scale))
gs = fig.add_gridspec(1,5)
dlr = 0.005
gs.update(left=dlr, right=1-dlr, top=0.975, bottom=0.35, wspace=0.15, hspace=0.0)
geo, prj = sfplt.getprojection(rotation=55-90, inclination=50) 
ax1 = fig.add_subplot(gs[0,0], projection=prj)
ax2 = fig.add_subplot(gs[0,1], projection=prj)
ax3 = fig.add_subplot(gs[0,2], projection=prj)
ax4 = fig.add_subplot(gs[0,3], projection=prj)
ax5 = fig.add_subplot(gs[0,4], projection=prj)
for ax in (ax1,ax2,ax3,ax4,ax5): ax3.set_global()

plm0, qlm0 = np.zeros(5), np.zeros(3)

sr = np.sqrt(2)
cij = np.exp(1j*np.deg2rad(0)) / sr
ccij = np.conjugate(cij)

cbaspect = 9
kwargs = dict(lvls=np.arange(0, 0.5+0.01, 0.1), cbtickintvl=2, cmap='RdPu', arrscale=4, cbaspect=cbaspect)
sfplt.plot_VSH(plm0, np.array([0,np.sqrt(2)*cij,0]), ax1, cblabel=r'$\norm*{\vb{Q}_0^0}$', **kwargs)
sfplt.plot_VSH(plm0, np.array([-ccij,0,cij]),        ax2, cblabel=r'$\norm*{\frac{1}{\sqrt{2}}\vb{Q}_0^0 + \text{c.c.}}$', **kwargs)

kwargs = dict(lvls=np.arange(0, 1.0+0.01, 0.25), cbtickintvl=2, cmap='GnBu', arrscale=7.5, cbaspect=cbaspect)
sfplt.plot_VSH(np.array([0,0,np.sqrt(2)*cij,0,0]), qlm0, ax3, cblabel=r'$\norm*{\vb{P}_0^0}$', **kwargs)
sfplt.plot_VSH(np.array([0,-ccij,0,cij,0]), qlm0,        ax4, cblabel=r'$\norm*{\frac{1}{\sqrt{2}}\vb{P}_2^1 + \text{c.c.}}$', **kwargs)
sfplt.plot_VSH(np.array([ccij,0,0,0,cij]), qlm0,         ax5, cblabel=r'$\norm*{\frac{1}{\sqrt{2}}\vb{P}_2^2 + \text{c.c.}}$', **kwargs)

for ax in (ax1,ax2,ax3,ax4,ax5):
    sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', color='k', fontsize=FS+1)

plt.savefig('VSH-modes.pdf', transparent=False, dpi=250)

