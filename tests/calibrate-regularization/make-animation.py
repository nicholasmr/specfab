# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
#from numpy import linalg as LA
#from scipy.optimize import minimize

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc

import warnings
warnings.filterwarnings("ignore")

#---------------------
# Setup
#---------------------

L = 20
Nt = 100 # Number of "integration time steps" (we consider the analytical solution here).

### Mode of deformation: unconfined vertical compression
epszz = -0.97 # target strain_zz
D = np.diag([0.5,0.5,-1]);
Dzz = D[2,2]; # size of compressive strain-rate
W = 0*D # Assume strain-rate tensor is in eigen frame

### Numerics
te = 1/Dzz # characteristic time, t_e
t_epszz = te*np.log(epszz + 1)
dt = t_epszz/Nt # time step size for thresshold strain epszz in "Nt" time steps 
tsteps = np.arange(Nt) # list of steps
epszz_t = np.exp(tsteps*dt/te)-1 # list of strains
print('Nt=%i, dt=%.4e s'%(Nt,dt))

#---------------------
# Fabric evolution
#---------------------

lm, nlm_len = sf.init(L)
nlm = np.zeros((2, nlm_len, Nt), dtype=np.complex128)
nlm[:,0,:] = 1/np.sqrt(4*np.pi) # Init with isotropy
nlm0 = nlm[0,:,0].copy()

M_LROT = sf.M_LROT(nlm0, D, W, 1, 0) # strain-rate assumed constant for calibration experiments
M_REG  = sf.M_REG(nlm0, D)

for ii in [0,1]:

    if ii == 0: M = M_LROT + M_REG 
    else:       M = M_LROT

    for tt in np.arange(1,Nt):

        nlm_prev = nlm[ii,:,tt-1].copy()
        nlm[ii,:,tt] = nlm_prev + dt*np.matmul(M, nlm_prev)

#---------------------
# Power spectrum, S(l)
#---------------------

nlm_dirac = nlm0.copy()
for ii, (l,m) in enumerate(lm.T): nlm_dirac[ii] = sp.sph_harm(m,l, 0,0)
Lrange = np.arange(0,L+1,2) # 0, 2, 4, 6, ...
Sl_dirac = np.array([sf.Sl(nlm_dirac, l) for l in Lrange])
Sl_dirac /= Sl_dirac[0] # normalize
print('S_dirac(l) (for l=0,2,...,L) = ', Sl_dirac)

#---------------------
# Plot results
#---------------------

# Plot frames
if 1:
    for tt in np.arange(Nt):
        
        ### Setup figure

        scale=1.0
        fig = plt.figure(figsize=(7*scale,3*scale))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.5,1,1], wspace=0.25, bottom=0.17, top=0.95, left=0.1, right=0.99)
        
        ### Plot power spectra

        ax = fig.add_subplot(gs[0,0])
        Sl_model_reg   = np.array([sf.Sl(nlm[0,:,tt], l) for l in Lrange]) 
        Sl_model_noreg = np.array([sf.Sl(nlm[1,:,tt], l) for l in Lrange]) 
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
        
        inclination = 50 # view angle
        rot0 = -90
        rot = -40 + rot0 
        prj = ccrs.Orthographic(rot, 90-inclination)
        geo = ccrs.Geodetic()     
        ax_ODF = [fig.add_subplot(gs[0, 1+ii], projection=prj) for ii in [0,1]]

        lvls = np.linspace(0,1,9)
        tickintvl = 4

        plot_ODF(nlm[0,:,tt], lm, ax=ax_ODF[0], cmap='Greys', cblabel=r'$n/N$ (ODF)', latres=40, lvls=lvls, tickintvl=tickintvl)
        ax_ODF[0].set_global()
        ax_ODF[0].set_title(r'$\bf{M}=\bf{M}_\mathrm{LROT} + \bf{M}_\mathrm{REG}$', fontsize=FS+1, pad=10)
        
        plot_ODF(nlm[1,:,tt], lm, ax=ax_ODF[1], cmap='Greys', cblabel=r'$n/N$ (ODF)', latres=40, lvls=lvls, tickintvl=tickintvl)
        ax_ODF[1].set_global()
        ax_ODF[1].set_title(r'$\bf{M}=\bf{M}_\mathrm{LROT}$', fontsize=FS+1, pad=10)

        ### Save figure

        fout = 'frames/calibrate-regularization-%03i.png'%(tt)
        print('Saving %s'%(fout))
        plt.savefig(fout, dpi=200)
    
# Make GIF
if 1:
    os.system(r'rm animation.gif')  
    os.system(r'ffmpeg -y -f image2 -framerate 40 -stream_loop 0 -i frames/calibrate-regularization-%03d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p animation.avi')
    os.system(r'ffmpeg -i animation.avi -vf "fps=15,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 animation.gif')
    os.system(r'rm animation.avi')  
    
