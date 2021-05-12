# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

import numpy as np
import scipy.special as sp
from netCDF4 import Dataset
import sys, os, copy, code # code.interact(local=locals())

#------------------
# Input arguments
#------------------

if len(sys.argv) != 3: 
    print('usage: %s nprime ugrad'%(sys.argv[0]))
    sys.exit(0)

exprref = sys.argv[2]
nprime = int(sys.argv[1])

# Options
latres = 40 # latitude resolution on S^2
inclination = 50 # view angle
rot0 = -90
rot = 1.4*rot0 # view angle

#------------------
# Load solution
#------------------

fh = Dataset('solutions/solution_n%i_%s.nc'%(nprime, exprref), mode='r')
loadvar = lambda field: np.array(fh.variables[field][:])

# Model config
Nt, dt = fh.getncattr('tsteps'), fh.getncattr('dt')

# Fabric state
lm, c   = loadvar('lm'), loadvar('c_re') + 1j*loadvar('c_im') 
a2_true = loadvar('a2_true')
a2      = loadvar('a2')

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

    dpi, scale = 200, 2.6
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

    lvls = np.arange(0,0.2+1e-5,0.025) # Contour lvls

    F = np.real(np.sum([ c[tt,ii]*sp.sph_harm(m, l, phi,theta) for ii,(l,m) in enumerate(lm) ], axis=0))
    F = np.divide(F, c[0,0]*np.sqrt(4*np.pi))
    hdistr = axdistr.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='max', cmap='Greys')
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = axdistr.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    cb1 = plt.colorbar(hdistr, ax=axdistr, fraction=0.04, aspect=19,  orientation='horizontal', pad=0.1, ticks=lvls[::2])   
    cb1.set_label(r'$n(\theta,\phi)/N(t=0)$', fontsize=FS)
    #
    ax = axdistr
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
        eigenValues, eigenVectors = np.linalg.eig(A)
#        idx = np.argsort(eigenValues)
#        eigenValues = eigenValues[idx,:]
        return (eigenValues)

    eigvals_true = eigen(a2_true)
    eigvals      = eigen(a2)

    steps = np.arange(len(eigvals[:,0]))

    lw = 1.5
    lwtrue = lw+0.5
    
    axeigvals.plot([tt,tt],[0,1],':k',lw=lw+1)
    #
    axeigvals.plot(steps,eigvals[:,0],'--r',label='$a_{1}$ (tens.)',lw=lw+0.5)
    axeigvals.plot(steps,eigvals[:,1],'--g',label='$a_{2}$ (tens.)',lw=lw+1.0)
    axeigvals.plot(steps,eigvals[:,2],'--b',label='$a_{3}$ (tens.)',lw=lw)
    #
    axeigvals.plot(steps,eigvals_true[:,0],'-r',label='$a_{1}$ (spec.)',lw=lwtrue+0.5)
    axeigvals.plot(steps,eigvals_true[:,1],'-g',label='$a_{2}$ (spec.)',lw=lwtrue+1.0)
    axeigvals.plot(steps,eigvals_true[:,2],'-b',label='$a_{3}$ (spec.)',lw=lwtrue)
    
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
    fout = 'solutions/solution_n%i_%s__%i.png'%(nprime,exprref, tt+1)
    print('Saving %s'%(fout))
    plt.savefig(fout,dpi=dpi)

