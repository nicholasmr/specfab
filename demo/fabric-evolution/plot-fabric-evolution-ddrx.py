# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020-2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp
from netCDF4 import Dataset

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import common as sfcom
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

import warnings
warnings.filterwarnings("ignore")

#------------------
# Input arguments
#------------------

if len(sys.argv) != 2: 
    print('usage: %s ugrad'%(sys.argv[0]))
    sys.exit(0)

exprref = sys.argv[1]

#------------------
# Load solution
#------------------

fh = Dataset('solutions/DDRX_%s.nc'%(exprref), mode='r')
loadvar = lambda field: np.array(fh.variables[field][:])

# Model config
Nt, dt = fh.getncattr('tsteps'), fh.getncattr('dt')

# Fabric state
lm, c   = loadvar('lm'), loadvar('c_re') + 1j*loadvar('c_im') 
lm = np.array(lm).T
a2_spec = loadvar('a2_true')
a2_tens = loadvar('a2')

#------------------
# Plot
#------------------

tsteps = [0, int(Nt*1/2), Nt-1] # time steps to plot
#tsteps = [0, Nt-1] # time steps to plot

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

rotation=55 
inclination=45
geo, prj = sfplt.getprojection(rotation=rotation, inclination=inclination)

for tt in tsteps:

    #----------------------
    # Figure setup
    #----------------------

    dpi, scale = 200, 2.8
    fig = plt.figure(figsize=(2.5*scale,1.3*scale))
    gs = gridspec.GridSpec(1,2, width_ratios=[1,1.2])
    a = 0.04
    gs.update(left=a, right=1-a, top=0.95, bottom=0.16, wspace=0.015*18, hspace=0.25)

    ax_ODF = plt.subplot(gs[0, 0], projection=prj)
    ax_ai  = plt.subplot(gs[0, 1])

    ax_ODF.set_global() # show entire S^2
                
    #----------------------
    # ODF
    #----------------------

    sfplt.plotODF(c[tt,:], lm, ax_ODF, cmap='Greys')
    sfplt.plotcoordaxes(ax_ODF, geo, axislabels='vuxi')

    # a4 eigentensors
    Wi = [sfcom.eigenframe4(c[tt,:], i=i, modelplane='xy')[0] for i in range(6)]
    code.interact(local=locals())
#    Q1,Q2,Q3,Q4,Q5,Q6, eigvals6 = sf.a4_eigentensors(c[tt,:])
#    _, W1 = np.linalg.eig(Q6)

    ms = 7
    kwargs = dict(ms=ms, markerfacecolor='none', markeredgewidth=1.0, transform=geo)
    ci = ['r','b','g','m','c','y']
    mrk = ['o','s','d','^','v','>']
    
    for ii in range(6):
        wi = Wi[ii]
        for kk in range(3):
            ei = wi[:,kk]
            sfplt.plotS2point(ax_ODF, +ei, marker=mrk[ii], markeredgecolor=ci[ii], **kwargs)
            sfplt.plotS2point(ax_ODF, -ei, marker=mrk[ii], markeredgecolor=ci[ii], **kwargs)

    #----------------------
    # Eigenvalues
    #----------------------
    
    def eigenvalues(A):
    
        N,_,_ = np.shape(A)
        ai = np.zeros((N,3))
        
        for nn in np.arange(N):
            ai_nn, _ = np.linalg.eig(A[nn,:,:])
            ai[nn,:] = ai_nn[np.flip(np.argsort(ai_nn))]
            
        return ai

    #----------

    lw = 1.5
    lwtrue = lw+0.2
    
    ax_ai.plot([tt,tt],[0,1],':k',lw=lw+1)

    ai_spec, ai_tens = eigenvalues(a2_spec), eigenvalues(a2_tens)
    steps = np.arange(len(ai_spec[:,0]))

    ax_ai.plot(steps, ai_spec[:,0], '-', color='#e31a1c', label='$a_{1}$ (spec.)', lw=lwtrue+0.25)
    ax_ai.plot(steps, ai_spec[:,1], '-', color='#33a02c', label='$a_{2}$ (spec.)', lw=lwtrue+0.50)
    ax_ai.plot(steps, ai_spec[:,2], '-', color='#1f78b4', label='$a_{3}$ (spec.)', lw=lwtrue)

    ax_ai.plot(steps, ai_tens[:,0], '--', color='#fb9a99', label='$a_{1}$ (tens.)', lw=lw+0.25)
    ax_ai.plot(steps, ai_tens[:,1], '--', color='#b2df8a', label='$a_{2}$ (tens.)', lw=lw+0.50)
    ax_ai.plot(steps, ai_tens[:,2], '--', color='#a6cee3', label='$a_{3}$ (tens.)', lw=lw)
    
    ax_ai.plot([0,1], [-tt,-tt], ':k', lw=lw)
    ax_ai.set_ylim([0,1])
    ax_ai.set_xlim([0,Nt])
    ax_ai.set_xlabel('time step')
    ax_ai.set_ylabel('$a_{i}$')
    ax_ai.grid()              
    ax_ai.legend(handlelength=1, ncol=2, labelspacing=0.3, fancybox=False, loc=2)
            
    #----------------------
    # Save figure
    #----------------------
    
    fout = 'solutions/DDRX_%s__%i.png'%(exprref, tt+1)
    print('Saving %s'%(fout))
    plt.savefig(fout,dpi=dpi)

