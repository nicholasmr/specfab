# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, rc

PLOT_FRAMES = 1
MAKE_GIF = 1

#---------------------
# Setup
#---------------------

L = 20 
Nt = 100 
strain_target = -0.97

mod = dict(type='ps', axis=2, T=1, r=0) # mode of deformation

#---------------------
# Fabric evolution
#---------------------

lm, nlm_len = sf.init(L)

# latrot only with calibrated high-L regularization
nlm_reg, *_ = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=+1, nu=1) 

# latrot only without regularization
nlm_noreg, *_ = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=+1, nu=None) 

#---------------------
# Plot results
#---------------------

# Plot frames
if PLOT_FRAMES:
    for tt in np.arange(Nt+1):
        
        ### Setup figure

        scale=0.9
        fig = plt.figure(figsize=(7*scale,3*scale))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.5,1,1], wspace=0.25, bottom=0.17, top=0.95, left=0.11, right=0.99)
        
        ### Plot power spectra

        Sl_dirac, Lrange, nlm_dirac = sfdsc.Sl_delta(lm[0,-1])

        ax = fig.add_subplot(gs[0,0])
        Sl_model_reg   = np.array([sf.Sl(nlm_reg[tt,:], l)   for l in Lrange]) 
        Sl_model_noreg = np.array([sf.Sl(nlm_noreg[tt,:], l) for l in Lrange]) 
        Sl_model_reg   /= Sl_model_reg[0] # normalize
        Sl_model_noreg /= Sl_model_noreg[0] # normalize
        h = ax.semilogy(Lrange, Sl_model_noreg, ls='-', c='#a50f15', label=r'$\bf{M}=\bf{M}_\mathrm{LROT}$')
        h = ax.semilogy(Lrange, Sl_model_reg,   ls='-', c='#006d2c', label=r'$\bf{M}=\bf{M}_\mathrm{LROT} + \bf{M}_\mathrm{REG}$')
            
        # Delta spectrum
        ax.semilogy(Lrange, Sl_dirac, '--', c='k', lw=1.5, label=r'$n(\bf{r})=\delta(\bf{r}-\bf{z})$')  
        
        # Set figure axes etc.
        ax.set_xlabel('$l$')
        ax.set_ylabel('$S(l)/S(0)$')
        ax.set_ylim([1e-3,10])
        ax.set_xticks(np.arange(0,21,2))
        ax.set_xlim([0,np.amax([10,L])])
        ax.grid()
        legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'handlelength':1.9, 'labelspacing':0.2}
        hleg = ax.legend(loc=3, fontsize=FS-0.5, **legkwargs)
        
        ### Plot ODFs
        
        geo, prj = sfplt.getprojection(rotation=45, inclination=50)
        ax_ODF = [fig.add_subplot(gs[0, 1+ii], projection=prj) for ii in [0,1]]

        lvlset = [np.linspace(0,1,9), lambda x,p:'%.1f'%x]

        sfplt.plotODF(nlm_reg[tt,:], lm, ax_ODF[0], lvlset=lvlset)
        sfplt.plotcoordaxes(ax_ODF[0], geo, axislabels='vuxi')
        ax_ODF[0].set_global()
        ax_ODF[0].set_title(r'$\bf{M}=\bf{M}_\mathrm{LROT} + \bf{M}_\mathrm{REG}$', fontsize=FS+1, pad=10)
        
        sfplt.plotODF(nlm_noreg[tt,:], lm, ax=ax_ODF[1], lvlset=lvlset)
        sfplt.plotcoordaxes(ax_ODF[1], geo, axislabels='vuxi')
        ax_ODF[1].set_global()
        ax_ODF[1].set_title(r'$\bf{M}=\bf{M}_\mathrm{LROT}$', fontsize=FS+1, pad=10)

        ### Save figure

        fout = 'frames/calibrate-regularization-%03i.png'%(tt)
        print('Saving %s'%(fout))
        plt.savefig(fout, dpi=200)
        plt.close()
    

if MAKE_GIF:
    os.system(r'rm regularization.gif')  
    os.system(r'ffmpeg -y -f image2 -framerate 40 -stream_loop 0 -i frames/calibrate-regularization-%03d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p regularization.avi')
    os.system(r'ffmpeg -i regularization.avi -vf "fps=15,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 regularization.gif')
    os.system(r'rm regularization.avi')  
    
