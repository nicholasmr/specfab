# N. M. Rathmann <rathmann@nbi.ku.dk>, 2024

"""
Plot effect of regularization on power spectrum
"""

import sys, os, copy, code # code.interact(local=locals())

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, rc

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
FS = 12
FS = sfplt.setfont_tex(fontsize=FS)
FSAX = FS+0

#---------------------
# Setup
#---------------------

Nt = 100 
strain_target = -0.97
mod = dict(type='ps', axis=2, T=1, r=0) # mode of deformation

#---------------------
# Fabric evolution
#---------------------

args   = (sf, mod, strain_target)
kwargs = dict(Nt=Nt, iota=+1)

L = 16
lm1, nlm_len = sf.init(L)
nlm1, F, *_       = sfint.lagrangianparcel(*args, **kwargs, nu=1) 
nlm1_noreg, *_ = sfint.lagrangianparcel(*args, **kwargs, nu=None) 

Ll = 8
lm2, nlm_len = sf.init(Ll)
nlm2, *_       = sfint.lagrangianparcel(*args, **kwargs, nu=1) 
nlm2_noreg, *_ = sfint.lagrangianparcel(*args, **kwargs, nu=None) 

Sl1_dirac, Lrange1, nlm1_dirac = sfdsc.Sl_delta(lm1[0,-1])
Sl2_dirac, Lrange2, nlm2_dirac = sfdsc.Sl_delta(lm2[0,-1])

for kk, tt in enumerate((20,60)):
#for kk, tt in enumerate((20,)):

    ### Setup figure

    scale=0.9
    fig = plt.figure(figsize=(4.5*scale,3.8*scale))
    gs = gridspec.GridSpec(2, 3, width_ratios=[3,1,1], wspace=0.12, hspace=1.0, \
        bottom=0.28, top=0.94, left=0.155, right=0.99)

    ax = fig.add_subplot(gs[:,0])
    
    geo, prj = sfplt.getprojection(rotation=60-2*90, inclination=50)
    ax_ODF = [fig.add_subplot(gs[int(np.floor(ii/2)), 1+ii%2], projection=prj) for ii in range(4)]
    for ii in range(4): ax_ODF[ii].set_global()

    ### Plot power spectra

    Sl1       = np.array([sf.Sl(nlm1[tt,:], l)       for l in Lrange1]) 
    Sl1_noreg = np.array([sf.Sl(nlm1_noreg[tt,:], l) for l in Lrange1])
    Sl2       = np.array([sf.Sl(nlm2[tt,:], l)       for l in Lrange2]) 
    Sl2_noreg = np.array([sf.Sl(nlm2_noreg[tt,:], l) for l in Lrange2]) 
     
    Sl1       /= Sl1[0]       # normalize
    Sl1_noreg /= Sl1_noreg[0] 
    Sl2       /= Sl2[0]       
    Sl2_noreg /= Sl2_noreg[0]

    h = ax.semilogy(Lrange1, Sl1,       ls='-', c=sfplt.c_dpurple, label=r'$L=%i$ + reg.'%(L))
    h = ax.semilogy(Lrange2, Sl2,       ls='-', c=sfplt.c_dgreen, label=r'$L=%i$ + reg.'%(Ll))

    ax.semilogy(Lrange1, Sl1_dirac, '--', c='k', lw=1.5, label=r'$n(\vu{r})=\delta(\vu{r}-\vu{z})$')

    h = ax.semilogy(Lrange1, Sl1_noreg, ls=':', c=sfplt.c_purple, label=r'$L=%i$'%(L))
    h = ax.semilogy(Lrange2, Sl2_noreg, ls=':', c=sfplt.c_green, label=r'$L=%i$'%(Ll))

    # Set figure axes etc.
    ax.set_xlabel('$l$')
    ax.set_ylabel('$S(l)/S(0)$')
    ax.set_ylim([1e-3,2])
    ax.set_xticks(np.arange(0,21,4))
    ax.set_xlim([0,np.amax([10,L])])
    ax.grid()
    legkwargs = {'frameon':False, 'fancybox':False, 'edgecolor':'k', 'columnspacing':1, 'handlelength':1.2, 'labelspacing':0.2}
    hleg = ax.legend(loc=3, fontsize=FS, ncol=2, bbox_to_anchor=(-0.25,-0.49), **legkwargs)

    ### Plot ODFs

    lvlset = [np.linspace(0,1,9), lambda x,p:'%.0f'%x]
    kw1 = dict(lvlset=lvlset, cblabel=r'$n(\vu{r})$', cblabelpad=0, cbtickintvl=8)

    ax = ax_ODF[0]
    sfplt.plotODF(nlm1[tt,:], lm1, ax=ax, cmap='Purples', **kw1)
    sfplt.plotcoordaxes(ax, geo, fontsize=FSAX, color='k', axislabels='vuxi')
    ax.set_title(r'$L=%i$ + reg.'%(L), fontsize=FS, pad=5)
    
    ax = ax_ODF[1]
    sfplt.plotODF(nlm1_noreg[tt,:], lm1, ax=ax, cmap='Purples', **kw1)
    sfplt.plotcoordaxes(ax, geo, fontsize=FSAX, color='k', axislabels='vuxi')
    ax.set_title(r'$L=%i$'%(L), fontsize=FS, pad=5)

    ax = ax_ODF[2]
    sfplt.plotODF(nlm2[tt,:], lm2, ax=ax, cmap='Greens', **kw1)
    sfplt.plotcoordaxes(ax, geo, fontsize=FSAX, color='k', axislabels='vuxi')
    ax.set_title(r'$L=%i$ + reg.'%(Ll), fontsize=FS, pad=5)
    
    ax = ax_ODF[3]
    sfplt.plotODF(nlm2_noreg[tt,:], lm2, ax=ax, cmap='Greens', **kw1)
    sfplt.plotcoordaxes(ax, geo, fontsize=FSAX, color='k', axislabels='vuxi')
    ax.set_title(r'$L=%i$'%(Ll), fontsize=FS, pad=5)

    ### Parcel plot
        
    W = 0.18
    ax1 = plt.axes([0.4 - kk*0.09, 0.29, W,W], projection='3d')
#    ax1.patch.set_visible(False)
    kwargs = dict(azim=35, axscale=1, axislabels=True, drawinit=True, fonttex=True, lw=0.75, colorinit=sfplt.c_dred)
    sfplt.plotparcel(ax1, F[tt], **kwargs)
    eps = sf.F_to_strain(F[tt])
    ax1.set_title('$\epsilon_{zz}=%.1f$'%(eps[-1,-1]), fontsize=FS, pad=-2)

    ### Save figure

    fout = 'regularization-powerspectrum-%i.pdf'%(tt)
    print('Saving %s'%(fout))
    plt.savefig(fout, dpi=200)
    plt.close()
                
