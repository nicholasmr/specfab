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
L = int(sys.argv[1])

# Options
latres = 40 # latitude resolution on S^2
inclination = 50 # view angle
rot = -90*1.4 # view angle

#------------------
# Load solution
#------------------

fh = Dataset('solutions/solution_n%i_%s.nc'%(L, exprref), mode='r')
loadvar = lambda field: np.array(fh.variables[field][:])

# Grain rheology
nprime, Ecc, Eca = fh.getncattr('nprime'), fh.getncattr('Ecc'), fh.getncattr('Eca')

# Model config
Nt, dt, nu, L = fh.getncattr('tsteps'), fh.getncattr('dt'), fh.getncattr('nu'), fh.getncattr('L')

# Fabric state
lm, c = loadvar('lm'), loadvar('c_re') + 1j*loadvar('c_im') 
eigvals  = loadvar('eigvals')
Eeiej    = loadvar('Eeiej')
Epijqij  = loadvar('Epijqij')
e1,e2,e3    = loadvar('e1'), loadvar('e2'), loadvar('e3')
p23,p12,p13 = loadvar('p23'),loadvar('p12'),loadvar('p13') 
q23,q12,q13 = loadvar('q23'),loadvar('q12'),loadvar('q13') 

# Velocity gradient tensor 
ugrad = np.array(fh.getncattr('ugrad')).reshape((3, 3))
eps = (ugrad + np.transpose(ugrad))/2
omg = (ugrad - np.transpose(ugrad))/2

#------------------
# Plot
#------------------

plot_tsteps = [0, int(Nt*1/2), Nt-1] # time steps to plot

lw0, lw1, lw2 = 3, 2.2, 1 

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

    dpi, scale = 200, 3.3
    fig = plt.figure(figsize=(2.1*scale,2.3*scale))
    gs = gridspec.GridSpec(2,2, height_ratios=[1,1.2], width_ratios=[1,1.2])
    a = 0.09
    gs.update(left=a, right=1-a/3, top=0.98, bottom=0.22, wspace=0.015*18, hspace=0.25)

    prj = ccrs.Orthographic(rot, 90-inclination)
    geo = ccrs.Geodetic()
    axdistr   = plt.subplot(gs[0, 0], projection=prj)
    axei      = plt.subplot(gs[0, 1], projection='3d')
    axeigvals = plt.subplot(gs[1, 0])
    axEvw     = plt.subplot(gs[1, 1])

    axdistr.set_global() # show entire S^2
                
    #----------------------
    # n(theta,phi) on S^2
    #----------------------

    lvls = np.linspace(0,1,5+1) # Contour lvls

    F = np.real(np.sum([ c[tt,ii]*sp.sph_harm(m, l, phi,theta) for ii,(l,m) in enumerate(lm) ], axis=0))
    hdistr = axdistr.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='max', cmap='Greys')
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = axdistr.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    cb1 = plt.colorbar(hdistr, ax=axdistr, fraction=0.04, aspect=19,  orientation='horizontal', pad=0.1)   
    cb1.set_label(r'$n(\theta,\phi)$', fontsize=FS)

    #----------------------
    # Eigen values
    #----------------------

    steps = np.arange(len(eigvals[:,0]))

    lw = 2
    axeigvals.plot([tt,tt],[0,1],':k',lw=lw)
    axeigvals.plot(steps,eigvals[:,0],'-r',label='$a_{1}$',lw=lw0)
    axeigvals.plot(steps,eigvals[:,1],'-g',label='$a_{2}$',lw=lw1)
    axeigvals.plot(steps,eigvals[:,2],'-b',label='$a_{3}$',lw=lw2)
    axeigvals.plot([0,1],[-tt,-tt],':k',lw=lw) 
    axeigvals.set_ylim([0,1])
    axeigvals.set_xlim([0, Nt+1])
    axeigvals.set_xlabel('time step')
    axeigvals.set_ylabel('$a_{i}$')
    axeigvals.grid()              
    axeigvals.legend(handlelength=1, ncol=1, labelspacing=0.3, fancybox=False, loc=2)
            
    #----------------------
    # Principal frame
    #----------------------

    axei.view_init(elev=90-inclination, azim=rot)
    axei.set_xlabel('$x$'); axei.set_xlim([-1,1])
    axei.set_ylabel('$y$'); axei.set_ylim([-1,1])
    axei.set_zlabel('$z$'); axei.set_zlim([0,1])
    plot_vec(axei,e1[tt,:],r'$\vb{e}_1$','r')
    plot_vec(axei,e2[tt,:],r'$\vb{e}_2$','g')
    plot_vec(axei,e3[tt,:],r'$\vb{e}_3$','b')
    lwpq=1
    plot_vec(axei,p13[tt,:],r'$\vb{p}_{13}$','y', lw=lwpq)
    plot_vec(axei,q13[tt,:],r'$\vb{q}_{13}$','y', ls='--', lw=lwpq)
    axei.legend(handlelength=1, loc=1, bbox_to_anchor=(1.17,1), fancybox=False)

    #----------------------
    # Enhancement factors
    #----------------------
    
    lblmk    = lambda ii,jj: '$E_{e_%i e_%i}$'%(ii+1,jj+1)
    lblmk_pq = lambda ii,jj: '$E_{p_{%i%i} q_{%i%i}}$'%(ii+1,jj+1, ii+1,jj+1)

    axEvw.semilogy(steps, Eeiej[:,0,0], 'r-',  label=lblmk(0,0), lw=lw0)
    axEvw.semilogy(steps, Eeiej[:,0,1], 'r--', label=lblmk(0,1), lw=lw0)
    axEvw.semilogy(steps, Eeiej[:,0,2], 'r:',  label=lblmk(0,2), lw=lw0)

    axEvw.semilogy(steps, Eeiej[:,1,0], 'g--', label=lblmk(1,0), lw=lw1)
    axEvw.semilogy(steps, Eeiej[:,1,1], 'g-',  label=lblmk(1,1), lw=lw1)
    axEvw.semilogy(steps, Eeiej[:,1,2], 'g-.', label=lblmk(1,2), lw=lw1)
        
    axEvw.semilogy(steps, Eeiej[:,2,0], 'b:',  label=lblmk(2,0), lw=lw2)
    axEvw.semilogy(steps, Eeiej[:,2,1], 'b-.', label=lblmk(2,1), lw=lw2)
    axEvw.semilogy(steps, Eeiej[:,2,2], 'b-',  label=lblmk(2,2), lw=lw2)    

    axEvw.semilogy(steps, Epijqij[:,2],  'y-', label=lblmk_pq(0,2), lw=lw1) 
    Eratio = np.divide(Eeiej[:,0,2],Epijqij[:,2])
    axEvw.semilogy(steps, Eratio, 'm-', label=lblmk(0,2)+'/'+lblmk_pq(0,2), lw=lw1)

    xlims=[0, Nt+1]
    lw=1
    axEvw.semilogy(xlims,[2.5,2.5],'--k',lw=lw) 
    axEvw.semilogy(xlims,[4.375,4.375],'--k',lw=lw) 
    axEvw.legend(loc=4, ncol=4, handlelength=2, columnspacing=1.2, labelspacing=0.3, fancybox=False, bbox_to_anchor=(1.05,-0.55))
    axEvw.set_xlim(xlims)
    axEvw.set_ylim([np.amin([1e-1, np.amin(Eeiej[:]), np.amin(Epijqij[:])]), np.amax([1e+1, np.amax(Eeiej[:]), np.amax(Epijqij[:]), np.amax(Eratio)])])
    axEvw.set_xlabel('time step')
    axEvw.set_ylabel('$E_{vw}$')
    axEvw.grid()
    
    #----------------------
    # Model config
    #----------------------
    
    props = dict(boxstyle='square', facecolor='wheat', alpha=0.5)
    textstr1 = '\n'.join(( r'"%s"'%(exprref.replace('_', '\_')), r"$n' = %i$"%(nprime), r'$E_{cc} = %.2f$'%(Ecc), r'$E_{ca} = %i$'%(Eca),  r'$L = %i$'%(L), r'$\nu = %.2e$'%(nu) ))
    axeigvals.text(-0.15, -0.54, textstr1, transform=axeigvals.transAxes, fontsize=FS, bbox=props)

    #----------------------
    # Save figure
    #----------------------
    fout = 'solutions/solution_n%i_%s__%i.png'%(nprime,exprref, tt+1)
    print('Saving %s'%(fout))
    plt.savefig(fout,dpi=dpi)

