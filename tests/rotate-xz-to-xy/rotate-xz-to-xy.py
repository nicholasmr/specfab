#!/usr/bin/python3
# Nicholas Rathmann

"""
Verify fortran routine rotate_nlm_xz2xy()
"""

import numpy as np
from scipy.spatial.transform import Rotation as R
from specfabpy import specfab as sf
from specfabpy import plotting as sfplt

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

L = 8
lm, nlm_len = sf.init(L)

FS = sfplt.setfont_tex(fontsize=11)

geoxz, prjxz = sfplt.getprojection(rotation=-90, inclination=90)
geoxy, prjxy = sfplt.getprojection(rotation=-90, inclination=0)

def plot(a2, fname):

    nlm = sf.a2_to_nlm(a2)
    
    fig = plt.figure(figsize=(5,2))
    gs = gridspec.GridSpec(1,3)
    ax1 = plt.subplot(gs[0,0], projection=prjxz)
    ax2 = plt.subplot(gs[0,1], projection=prjxz)
    ax3 = plt.subplot(gs[0,2], projection=prjxy)
    axes = (ax1,ax2,ax3)
    
    for ax in axes:
        ax.set_global()
        
        
    kwargs = dict(lvlset=(np.linspace(0.0, 0.3, 8), lambda x,p:'%.1f'%x), showcb=False)
    nlm_rot = sf.rotate_nlm_xz2xy(nlm)
    sfplt.plotODF(nlm,     lm, ax1, **kwargs)
    sfplt.plotODF(nlm_rot, lm, ax2, **kwargs)
    sfplt.plotODF(nlm_rot, lm, ax3, **kwargs)

    ax1.set_title(r'\texttt{original (xz)}')
    ax2.set_title(r'\texttt{rotated (xz)}')
    ax3.set_title(r'\texttt{rotated (xy)}')

    kwargs = dict(axislabels='vuxi', negaxes=False, color='k', fontsize=FS+1.5)
    sfplt.plotcoordaxes(ax1, geoxz, **kwargs)
    sfplt.plotcoordaxes(ax2, geoxz, **kwargs)
    sfplt.plotcoordaxes(ax3, geoxy, **kwargs)
        
    plt.tight_layout()
    plt.savefig('%s.png'%(fname), dpi=150)
    
    
a2_gd = np.diag([0.0,0.5,0.5])
M = R.from_rotvec(np.deg2rad(20) * np.array([0, 1, 0])).as_matrix()
a2 = np.matmul(M.T, np.matmul(a2_gd, M))
plot(a2, 'girdle')

a2_sm = np.diag([0.0,0.2,0.8])
M = R.from_rotvec(np.deg2rad(20) * np.array([0, 1, 0])).as_matrix()
a2 = np.matmul(M.T, np.matmul(a2_sm, M))
plot(a2, 'smax')
