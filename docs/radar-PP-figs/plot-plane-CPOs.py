#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import copy, os, sys, code # code.interact(local=locals())
import numpy as np
sys.path.insert(0, '../../demo')
os.system('mkdir -p ./frames')

from specfabpy import specfabpy as sf
from header import *

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc, colors
import matplotlib.gridspec as gridspec
import scipy.special as sp

### Plot

inclination = 45 # view angle
rot0 = 1 * -90 
rot = rot0 - 35 # view angle

prj = ccrs.Orthographic(rot, 90-inclination)
geo = ccrs.Geodetic()

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

def set_axis_labels(axlist, c='#99000d'):
    FSAX = FS+1
    for ax in axlist:
        ax.text(rot0-40, 88, r'$\vb{z}$', color=c, horizontalalignment='left', transform=geo, fontsize=FSAX)
        ax.text(rot0-96, -12, r'$\vb{m}_1$', color=c, horizontalalignment='left', transform=geo, fontsize=FSAX)
        ax.text(rot0-6, -5, r'$\vb{m}_2$', color=c, horizontalalignment='left', transform=geo, fontsize=FSAX)

def a2planar(lam1, dlam): 
    return np.diag([lam1, lam1+dlam, 1-dlam-2*lam1])

axlist, fig, fraction, aspect = setup_fig()

lvls = np.linspace(0.0,0.4,9)
tickintvl = 4

lm, nlm_len = sf.init(2)

for ii, dlam in enumerate([0, 0.25, 0.50, 0.75, 1]):
    plot_ODF(sf.a2_to_nlm(a2planar(0,dlam)),lm, ax=axlist[ii], cmap='Greys', cblabel='$n/N$ (ODF)', lvls=lvls, tickintvl=tickintvl, aspect=aspect, fraction=fraction)
    axlist[ii].set_title(r'$\Delta\lambdaup=%.2f$'%(dlam), fontsize=FS)

set_axis_labels(axlist)

plt.savefig('plan-CPOs.png', transparent=True, dpi=150)
