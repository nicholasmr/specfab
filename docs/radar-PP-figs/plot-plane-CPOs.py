#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import copy, os, sys, code # code.interact(local=locals())
import numpy as np
os.system('mkdir -p ./frames')

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc, colors
import matplotlib.gridspec as gridspec
import scipy.special as sp

### Plot

geo, prj = sfplt.getprojection(rotation=55, inclination=45)

def setup_fig():

    scale = 2.0
    fig = plt.figure(figsize=(5/2*1.2*scale,0.93*scale))
    gs = fig.add_gridspec(1,5)
    al = 0.025
    ar = 0.02
    gs.update(left=al, right=1-ar, top=0.88, bottom=0.25, wspace=0.3, hspace=0.4)

    ax1 = fig.add_subplot(gs[0,0], projection=prj); ax1.set_global(); 
    ax2 = fig.add_subplot(gs[0,1], projection=prj); ax2.set_global(); 
    ax3 = fig.add_subplot(gs[0,2], projection=prj); ax3.set_global(); 
    ax4 = fig.add_subplot(gs[0,3], projection=prj); ax4.set_global(); 
    ax5 = fig.add_subplot(gs[0,4], projection=prj); ax5.set_global(); 
    axlist = [ax1,ax2,ax3,ax4,ax5]

    fraction = 0.08
    aspect = 8

    return axlist, fig, fraction, aspect

def a2planar(lam1, dlam): 
    return np.diag([lam1, lam1+dlam, 1-dlam-2*lam1])

axlist, fig, fraction, aspect = setup_fig()

lm, nlm_len = sf.init(2)
#lvlset = 'iso-up'
lvlset = [np.linspace(0.0,0.4,9), lambda x,p:'%.1f'%x]

for ii, dlam in enumerate([0, 0.25, 0.50, 0.75, 1]):
    sfplt.plotODF(sf.a2_to_nlm(a2planar(0,dlam)), lm, axlist[ii], cblabel='MODF', lvlset=lvlset, cbaspect=aspect, cbfraction=fraction)
    sfplt.plotcoordaxes(axlist[ii], geo, axislabels=[r'$\vb{m}_1$',r'$\vb{m}_2$',r'$\vb{z}$'])
    axlist[ii].set_title(r'$\Delta\lambdaup=%.2f$'%(dlam), fontsize=FS)

plt.savefig('plane-CPOs.png', transparent=True, dpi=150)
