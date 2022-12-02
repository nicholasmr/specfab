# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

import numpy as np
import scipy.special as sp
import sys, os, copy, code # code.interact(local=locals())

sys.path.insert(0, '..')
from header import *
from specfabpy import specfabpy as sf

import warnings
warnings.filterwarnings("ignore")

### Options

latres = 40 # latitude resolution on S^2
#latres = 20 # latitude resolution on S^2
inclination = 50 # view angle
rot0 = -90
rot = 1.4*rot0 # view angle

### Numerics

Nt = 100 # Number of time steps
dt = 0.02 # Time-step size (gives a vertical strain of -0.98 for experiment "uc_zz")
L  = 8

### Process parameters

Gamma0 = 2e1 
Lam = 2e-1

### Integrate

lm, nlm_len = sf.init(L) 
nlm = np.zeros((3,3,Nt,nlm_len), dtype=np.complex64) # Array of expansion coefficients
nlm[:,:,0,0] = 1/np.sqrt(4*np.pi) # Normalized ODF at t=0

Mtypes = ['LROT','DDRX','CDRX']
#Mtypes = ['DDRX',]

for ii, Mtype in enumerate(Mtypes):

    if Mtype == 'CDRX':
        nlm[ii,:,0,:] = nlm[0,:,-1,:] # use LROT end state as starting point for CDRX

    for jj, expr in enumerate(['uc_zz', 'cc_zx', 'ss_xz']):
        
        ## Velocity gradient
        if expr == 'uc_zz': ugrad = np.diag([.5, .5, -1])
        if expr == 'cc_zx': ugrad = np.diag([1, 0, -1])
        if expr == 'ss_xz': ugrad = 2*np.array([[0,0,1], [0,0,0], [0,0,0]])

        D = (ugrad+np.transpose(ugrad))/2 # Symmetric part (strain-rate)
        W = (ugrad-np.transpose(ugrad))/2 # Anti-symmetric part (spin)
        S = D.copy() # stress tensor for DDRX
        
        ### Euler integration
        for tt in np.arange(1,Nt):
            nlm_prev = nlm[ii,jj,tt-1,:]
            if Mtype == 'LROT': M = sf.M_LROT(nlm_prev, D, W, 1, 0) + sf.M_REG(nlm_prev, D)
            if Mtype == 'DDRX': M = Gamma0 * sf.M_DDRX(nlm_prev, S)
            if Mtype == 'CDRX': M = Lam*sf.M_CDRX(nlm_prev)
            nlm[ii,jj,tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev)
        
### Plot

theta = np.linspace(0,   np.pi,   latres) # CO-LAT 
phi   = np.linspace(0, 2*np.pi, 2*latres) # LON
phi, theta = np.meshgrid(phi, theta) # gridded 
lon, colat = phi, theta
lat = np.pi/2-colat

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

FS = 13
rc('font',**{'family':'serif','sans-serif':['Times'],'size':FS})
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage{txfonts}'

def plot_vec(ax, v, lbl, color, ls='-', lw=2):
    ax.plot([0, +v[0]],[0, +v[1]],[0,+v[2]], color=color, ls=ls, lw=lw, label=lbl)
    ax.plot([0, -v[0]],[0, -v[1]],[0,-v[2]], color=color, ls=ls, lw=lw)

prj = ccrs.Orthographic(rot, 90-inclination)
geo = ccrs.Geodetic()

for ii, Mtype in enumerate(Mtypes):

    # Plot frames
    if 0: 
        for tt in np.arange(0,Nt):
        #for tt in np.arange(Nt-1,Nt):
        
            ### Figure setup
            scale = 2.6
            fig = plt.figure(figsize=(3/2*1.6*scale,1.00*scale))
            gs = fig.add_gridspec(1,3)
            al = 0.04
            ar = 0.02
            gs.update(left=al, right=1-ar, top=0.85, bottom=0.20, wspace=0.3, hspace=0.4)

            ax1 = fig.add_subplot(gs[0,0], projection=prj); ax1.set_global(); 
            ax2 = fig.add_subplot(gs[0,1], projection=prj); ax2.set_global(); 
            ax3 = fig.add_subplot(gs[0,2], projection=prj); ax3.set_global(); 
            axlist = [ax1,ax2,ax3]

            ### Plot
            lvls = np.linspace(0,1,9)
            tickintvl = 4
            plot_ODF(nlm[ii,0,tt,:], lm, ax=ax1, cmap='Greys', cblabel=r'$\psi/N$ (ODF)', lvls=lvls, tickintvl=tickintvl)
            plot_ODF(nlm[ii,1,tt,:], lm, ax=ax2, cmap='Greys', cblabel=r'$\psi/N$ (ODF)', lvls=lvls, tickintvl=tickintvl)
            plot_ODF(nlm[ii,2,tt,:], lm, ax=ax3, cmap='Greys', cblabel=r'$\psi/N$ (ODF)', lvls=lvls, tickintvl=tickintvl)
            ax1.set_title(r'%s Unconfined pure shear'%(r'{\Large \bf a}\,\,'), fontsize=FS, pad=10)
            ax2.set_title(r'%s Confined pure shear'%(r'{\Large \bf b}\,\,'), fontsize=FS, pad=10)
            ax3.set_title(r'%s Simple shear'%(r'{\Large \bf c}\,\,'), fontsize=FS, pad=10)
                
            ### Save
            fout = 'frames/fabric-evolution-%s-animation-%03i.png'%(Mtype,tt)
            print('Saving %s'%(fout))
            plt.savefig(fout, dpi=200)
            plt.close()

    # Make GIF
    if 1:
        os.system('rm animation-%s.gif'%(Mtype))  
        os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/fabric-evolution-%s-animation-%s.png -vcodec libx264 -crf 20  -pix_fmt yuv420p animation-%s.avi'%(Mtype,'%03d',Mtype))
        os.system('ffmpeg -i animation-%s.avi -vf "fps=20,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 animation-%s.gif'%(Mtype,Mtype))
        os.system('rm animation-%s.avi'%(Mtype))  

