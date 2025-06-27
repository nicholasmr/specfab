#!/usr/bin/python3
# Nicholas Rathmann

"""
Verify fortran routines for ideal nlm states
"""

import numpy as np
from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

lm, nlm_len = sf.init(20)
FS = sfplt.setfont_tex(fontsize=11)
geo, prj = sfplt.getprojection(rotation=-45, inclination=50)

def plot(m, L, fname):

    nlm_iso    = sf.nlm_isotropic(L)
    nlm_smax   = sf.nlm_singlemax(m,L)
    nlm_girdle = sf.nlm_girdle(m,L)

    fig = plt.figure(figsize=(4,2))
    gs = gridspec.GridSpec(1,3)
    ax1 = plt.subplot(gs[0,0], projection=prj)
    ax2 = plt.subplot(gs[0,1], projection=prj)
    ax3 = plt.subplot(gs[0,2], projection=prj)
#    axes = (ax1,ax2,ax3)        
        
    kwargs = dict(lvlset=(np.linspace(0.0, 0.3, 8), lambda x,p:'%.1f'%x), showcb=False)
    sfplt.plotODF(nlm_iso,    lm, ax1, **kwargs)
    sfplt.plotODF(nlm_smax,   lm, ax2, **kwargs)
    sfplt.plotODF(nlm_girdle, lm, ax3, **kwargs)

    for ax in (ax1,ax2,ax3): 
        ax.set_global()
        sfplt.plotcoordaxes(ax, geo)
        
    plt.tight_layout()
    plt.savefig('cpo-ideal-L%i-%s.png'%(L, fname), dpi=150, pad_inches=0.1, bbox_inches='tight')

for L in (2,4,8,12):
    plot(sfdsc.sph2cart(0,0), L, 'z')
    plot(sfdsc.sph2cart(np.pi/2,-np.pi/4), L, 'xy')

