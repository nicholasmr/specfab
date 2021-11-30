# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

import numpy as np
import scipy.special as sp
from netCDF4 import Dataset
import sys, os, copy, code # code.interact(local=locals())

sys.path.insert(0, '..')
from header import *

#------------------
# Input arguments
#------------------

if len(sys.argv) != 2: 
    print('usage: %s ugrad'%(sys.argv[0]))
    sys.exit(0)

exprref = sys.argv[1]

# Options
latres = 40 # latitude resolution on S^2
inclination = 50 # view angle
rot0 = -90
rot = 1.4*rot0 # view angle

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

plot_tsteps = [0, int(Nt*1/2), Nt-1] # time steps to plot
#plot_tsteps = [0, Nt-1] # time steps to plot

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

for tt in plot_tsteps:

    #----------------------
    # Figure setup
    #----------------------

    dpi, scale = 200, 2.8
    fig = plt.figure(figsize=(2.5*scale,1.3*scale))
    gs = gridspec.GridSpec(1,2, width_ratios=[1,1.2])
    a = 0.04
    gs.update(left=a, right=1-a, top=0.95, bottom=0.16, wspace=0.015*18, hspace=0.25)

    prj = ccrs.Orthographic(rot, 90-inclination)
    geo = ccrs.Geodetic()
    axdistr   = plt.subplot(gs[0, 0], projection=prj)
    axeigvals = plt.subplot(gs[0, 1])

    axdistr.set_global() # show entire S^2
                
    #----------------------
    # n(theta,phi) on S^2
    #----------------------

    ax = axdistr
    ODF = c[tt,:]/(np.sqrt(4*np.pi)*c[tt,0])
    plot_ODF(c[tt,:], lm, ax=ax, cmap='Greys', cblabel=r'$\psi/N$ (ODF)')
    ax.plot([0],[90],'k.', transform=geo) # z dir
    ax.plot([rot0],[0],'k.', transform=geo) # y dir
    ax.plot([rot0-90],[0],'k.', transform=geo) # x dir
    ax.text(rot0-80, 78, r'$\vu{z}$', horizontalalignment='left', transform=geo)
    ax.text(rot0-8, -8, r'$\vu{y}$', horizontalalignment='left', transform=geo)
    ax.text(rot0-90+5, -8, r'$\vu{x}$', horizontalalignment='left', transform=geo)

    #----------------------
    # Eigen values
    #----------------------
    
    def eigen(A):
        N,_,_ = np.shape(A)
        ai = np.zeros((N,3))
        for nn in np.arange(N):
            ai_nn, _ = np.linalg.eig(A[nn,:,:])
            ai[nn,:] = ai_nn[np.flip(np.argsort(ai_nn))]
        return ai

    ai_spec, ai_tens = eigen(a2_spec), eigen(a2_tens)

    steps = np.arange(len(ai_spec[:,0]))

    lw = 1.5
    lwtrue = lw+0.2
    
    axeigvals.plot([tt,tt],[0,1],':k',lw=lw+1)
    #
    axeigvals.plot(steps,ai_spec[:,0],'-', color='#e31a1c', label='$a_{1}$ (spec.)',lw=lwtrue+0.25)
    axeigvals.plot(steps,ai_spec[:,1],'-', color='#33a02c', label='$a_{2}$ (spec.)',lw=lwtrue+0.50)
    axeigvals.plot(steps,ai_spec[:,2],'-', color='#1f78b4', label='$a_{3}$ (spec.)',lw=lwtrue)
    #
    axeigvals.plot(steps,ai_tens[:,0],'--', color='#fb9a99', label='$a_{1}$ (tens.)',lw=lw+0.25)
    axeigvals.plot(steps,ai_tens[:,1],'--', color='#b2df8a', label='$a_{2}$ (tens.)',lw=lw+0.50)
    axeigvals.plot(steps,ai_tens[:,2],'--', color='#a6cee3', label='$a_{3}$ (tens.)',lw=lw)
    
    axeigvals.plot([0,1],[-tt,-tt],':k',lw=lw) 
    axeigvals.set_ylim([0,1])
    axeigvals.set_xlim([0, Nt+1])
    axeigvals.set_xlabel('time step')
    axeigvals.set_ylabel('$a_{i}$')
    axeigvals.grid()              
    axeigvals.legend(handlelength=1, ncol=2, labelspacing=0.3, fancybox=False, loc=2)
            
    #----------------------
    # Save figure
    #----------------------
    fout = 'solutions/DDRX_%s__%i.png'%(exprref, tt+1)
    print('Saving %s'%(fout))
    plt.savefig(fout,dpi=dpi)

