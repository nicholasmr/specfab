# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

"""
Supporting routines for CPO state space plotting etc.
These are not really for the end-user, but private Nicholas' tools.
"""

import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.colors

from . import discrete as sfdsc
from . import plotting as sfplt

c_iso    = sfplt.c_vdgray
c_planar = '#8c510a' 
c_unidir = '#01665e'
c_circle = '#762a83'

cl_planar = '#dfc27d'
cl_unidir = '#7fcdbb'
#cl_circle = 'tab:red'

cvl_unidir = matplotlib.colors.to_hex(np.array([199,234,229])/255)
cvl_planar = matplotlib.colors.to_hex(np.array([246,232,195])/255)
cvl_circle = '#fddbc7'

#-------------------------
# For state space plots
#-------------------------

def plot_nlm_cases(ax, fs, ms=7.5, norm=1/np.sqrt(4*np.pi), \
                    show_circle=True, isolbl_above=True, dy0=0.075, dx0_unidir=-0.03, dx0_planar=0.0,
                    c_planar=c_planar, c_circle=c_circle, c_unidir=c_unidir, c_iso=c_iso):

    nl0_unidir, nl0_planar, nl0_circle = sfdsc.nlm_ideal_cases(norm=norm)
    dy = dy0/norm
    
    kwargs_mrk = dict(marker='o', ms=ms, ls='none', zorder=20)
    kwargs_lbl = dict(ha='center', va='center', ma='center', fontsize=fs, zorder=20)
    
    ax.plot(0, 0, c=c_iso, **kwargs_mrk)
    plt.text(0, +dy if isolbl_above else -dy, r'{\bf Isotropic}', color=c_iso, **kwargs_lbl)

    ax.plot(*nl0_unidir, c=c_unidir, **kwargs_mrk)
    plt.text(nl0_unidir[0]+dx0_unidir/norm, nl0_unidir[1]+dy, '{\\bf Unidirectional}', color=c_unidir, **kwargs_lbl)

    ax.plot(*nl0_planar, c=c_planar, **kwargs_mrk)
    plt.text(nl0_planar[0]+dx0_planar/norm, nl0_planar[1]+dy, '{\\bf Planar}', color=c_planar, **kwargs_lbl)

    if show_circle:
        ax.plot(*nl0_circle, c=c_circle, **kwargs_mrk)
        plt.text(nl0_circle[0], nl0_circle[1]-dy, '{\\bf Circle}', color=c_circle, **kwargs_lbl)

def nlm_isvalid_grid(xlims, ylims, resx, resy):

    x = np.linspace(xlims[0], xlims[1], resx)
    y = np.linspace(ylims[0], ylims[1], resy)
    xv, yv = np.meshgrid(x, y, indexing='xy')
    
    isvalid = np.reshape(sfdsc.sf__.nlm_isvalid(xv.flatten(), yv.flatten()), (resy, resx))
    
    return isvalid, x, y
    
    
def Eij_statemap_ice(grain_params, xlims, ylims, resx=200, resy=200, ignore_isvalid=False):

    isvalid, x,y = nlm_isvalid_grid(xlims, ylims, resx, resy)

    Ezz = np.zeros((resy, resx)) 
    Exz = np.zeros((resy, resx)) 

    m, t = np.array([0,0,1]), np.array([1,0,0])
    mm, mt = np.tensordot(m,m, axes=0), np.tensordot(m,t, axes=0)
    tau_mm = np.identity(3) - 3*mm
    tau_mt = mt + mt.T

    for ii, y_ in enumerate(y):
        for jj, x_ in enumerate(x): 
        
            if ignore_isvalid or isvalid[ii,jj]:
                nlm_ = np.zeros((sfdsc.sf__.L4len), dtype=np.complex64)
                nlm_[0], nlm_[sfdsc.sf__.I20], nlm_[sfdsc.sf__.I40] = 1, x_, y_
                Ezz[ii,jj] = sfdsc.sf__.Evw_tranisotropic(nlm_, m,m,tau_mm, *grain_params)
                Exz[ii,jj] = sfdsc.sf__.Evw_tranisotropic(nlm_, m,t,tau_mt, *grain_params)
            else:
                Ezz[ii,jj] = np.nan
                Exz[ii,jj] = np.nan
    
    return (Ezz,Exz, x,y, isvalid)
    
