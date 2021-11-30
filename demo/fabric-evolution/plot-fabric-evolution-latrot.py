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

fh = Dataset('solutions/LATROT_%s.nc'%(exprref), mode='r')
loadvar = lambda field: np.array(fh.variables[field][:])

# Model config
Nt, dt, nu, L = fh.getncattr('tsteps'), fh.getncattr('dt'), fh.getncattr('nu'), fh.getncattr('L')

# Grain parameters
Eca_lin,  Ecc_lin  = fh.getncattr('Eca_lin'),  fh.getncattr('Ecc_lin')
Eca_nlin, Ecc_nlin = fh.getncattr('Eca_nlin'), fh.getncattr('Ecc_nlin')
alpha_lin, alpha_nlin = fh.getncattr('alpha_lin'), fh.getncattr('alpha_nlin')

# Fabric state
lm, c = loadvar('lm'), loadvar('c_re') + 1j*loadvar('c_im') 
lm = np.array(lm).T
eigvals = loadvar('eigvals')
e1,e2,e3 = loadvar('e1'), loadvar('e2'), loadvar('e3')
p1,p2,p3 = loadvar('p1'), loadvar('p2'), loadvar('p3')

Eeiej_lin, Eeiej_nlin = loadvar('Eeiej_lin'), loadvar('Eeiej_nlin') 
Epipj_lin, Epipj_nlin = loadvar('Epipj_lin'), loadvar('Epipj_nlin') 

# Velocity gradient tensor 
ugrad = np.array(fh.getncattr('ugrad')).reshape((3, 3))
eps = (ugrad + np.transpose(ugrad))/2
omg = (ugrad - np.transpose(ugrad))/2

#------------------
# Plot
#------------------

plot_tsteps = [0, int(Nt*1/2), Nt-1] # time steps to plot
#plot_tsteps = [Nt-1,] # time steps to plot

lw0, lw1, lw2 = 2.5,2.25,2.0

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
    fig = plt.figure(figsize=(3/2*2.1*scale,2.3*scale))
    gs = gridspec.GridSpec(2,3, height_ratios=[1,1.2], width_ratios=[1,1,1])
    a = 0.06
    gs.update(left=a, right=1-a/3, top=0.98, bottom=0.20, wspace=0.015*18, hspace=0.35)

    prj = ccrs.Orthographic(rot, 90-inclination)
    geo = ccrs.Geodetic()
    axdistr   = plt.subplot(gs[0, 0], projection=prj)
    axei      = plt.subplot(gs[0, 1], projection='3d')
    axeigvals = plt.subplot(gs[0, 2])
    axEvw_Mix = plt.subplot(gs[1, 1])
    axEvw_Sac = plt.subplot(gs[1, 2])

    axdistr.set_global() # show entire S^2
                
    #----------------------
    # n(theta,phi) on S^2
    #----------------------

    ax = axdistr
    ODF = c[tt,:]/(np.sqrt(4*np.pi)*c[tt,0])
    plot_ODF(c[tt,:], lm, ax=ax, cmap='Greys', cblabel=r'$\psi/N$ (ODF)')
    ax.plot([0],[90],'X', color=c_brown, transform=geo) # z dir
    ax.plot([rot0],[0],'X', color=c_brown, transform=geo) # y dir
    ax.plot([rot0-90],[0],'X', color=c_brown, transform=geo) # x dir
    ax.text(rot0-80, 78, r'$\vu{z}$', color=c_brown, horizontalalignment='left', transform=geo)
    ax.text(rot0-8, -8, r'$\vu{y}$', color=c_brown, horizontalalignment='left', transform=geo)
    ax.text(rot0-90+5, -8, r'$\vu{x}$', color=c_brown, horizontalalignment='left', transform=geo)

    #----------------------
    # Eigen values
    #----------------------

    steps = np.arange(len(eigvals[:,0]))

    lw = 2
    axeigvals.plot([tt,tt],[0,1],':k',lw=lw)
    axeigvals.plot(steps,eigvals[:,0],'-', c=c_red,   label='$a_{1}$',lw=lw0)
    axeigvals.plot(steps,eigvals[:,1],'-', c=c_green, label='$a_{2}$',lw=lw1)
    axeigvals.plot(steps,eigvals[:,2],'-', c=c_blue,  label='$a_{3}$',lw=lw2)
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
    plot_vec(axei,e1[tt,:],r'$\vb{e}_1$', c_red)
    plot_vec(axei,e2[tt,:],r'$\vb{e}_2$', c_green)
    plot_vec(axei,e3[tt,:],r'$\vb{e}_3$', c_blue)
    lwpq=1
    plot_vec(axei,p1[tt,:],r'$\vb{p}_{1}$', c_red,   ls='--', lw=lwpq)
    plot_vec(axei,p2[tt,:],r'$\vb{p}_{2}$', c_green, ls='--', lw=lwpq)
    plot_vec(axei,p3[tt,:],r'$\vb{p}_{3}$', c_blue,  ls='--', lw=lwpq)
    axei.legend(handlelength=1, loc=1, bbox_to_anchor=(1.17,1), fancybox=False)

    #----------------------
    # Enhancement factors
    #----------------------
    
    lblmk    = lambda ii,jj: '$E_{e_%i e_%i}$'%(ii+1,jj+1)
    lblmk_pq = lambda ii,jj: '$E_{p_{%i%i} q_{%i%i}}$'%(ii+1,jj+1, ii+1,jj+1)

    def plot_enh(ax, Eeiej, Epipj):
        
        ax.semilogy(steps, Eeiej[:,0,0], '-', c=c_red,  label=lblmk(0,0), lw=lw0)
        ax.semilogy(steps, Eeiej[:,1,1], '-', c=c_green, label=lblmk(1,1), lw=lw1)
        ax.semilogy(steps, Eeiej[:,2,2], '-', c=c_blue, label=lblmk(2,2), lw=lw2)    
        
        ax.semilogy(steps, Eeiej[:,0,1], ':', c=c_lred,   label=lblmk(0,1), lw=lw0)
        ax.semilogy(steps, Eeiej[:,1,2], ':', c=c_lgreen, label=lblmk(1,2), lw=lw2)
        ax.semilogy(steps, Eeiej[:,0,2], ':', c=c_lblue,  label=lblmk(0,2), lw=lw1)

# Don't plot, Eij is symmetric
#        ax.semilogy(steps, Eeiej[:,1,0], 'g--', label=lblmk(1,0), lw=lw1)
#        ax.semilogy(steps, Eeiej[:,2,0], 'b:',  label=lblmk(2,0), lw=lw2)
#        ax.semilogy(steps, Eeiej[:,2,1], 'b-.', label=lblmk(2,1), lw=lw2)

        ax.semilogy(steps, Epipj[:,0,1],  '-', c=c_gray, label=lblmk_pq(0,2), lw=lw2) 
        Eratio = np.divide(Eeiej[:,0,1],Epipj[:,0,1])
        ax.semilogy(steps, Eratio, '--', c=c_gray, label=lblmk(0,2)+'/'+lblmk_pq(0,2), lw=lw2)

        xlims=[0, Nt+1]
        lw=1
        ax.semilogy(xlims,[2.5,2.5],'--k',lw=lw) 
        ax.semilogy(xlims,[4.375,4.375],'--k',lw=lw) 
        ax.set_xlim(xlims)
        ax.set_ylim([np.amin([1e-1, np.amin(Eeiej[:]), np.amin(Epipj[:])]), np.amax([1e+1, np.amax(Eeiej[:]), np.amax(Epipj[:]), np.amax(Eratio)])])
        ax.set_xlabel('time step')
        ax.set_ylabel('$E_{vw}$')
        ax.grid()
    
    
    plot_enh(axEvw_Sac, Eeiej_nlin, Epipj_nlin)
    plot_enh(axEvw_Mix, Eeiej_lin,  Epipj_lin)
    
    axEvw_Sac.set_title(r'Nonlinear Sachs')
    axEvw_Mix.set_title(r'Linear mixed Taylor--Sachs')
    
    from matplotlib.offsetbox import AnchoredText
    axEvw_Sac.add_artist(AnchoredText('\n'.join( (r"$n' = %i$"%(3), r'$E_{cc} = %.1f$'%(Ecc_nlin), r'$E_{ca} = %.0e$'%(Eca_nlin), r'$\alpha = %.4f$'%(alpha_nlin)) ), loc='lower left', frameon=True,))
    axEvw_Mix.add_artist(AnchoredText('\n'.join( (r"$n' = %i$"%(1), r'$E_{cc} = %.1f$'%(Ecc_lin),  r'$E_{ca} = %.0e$'%(Eca_lin),  r'$\alpha = %.4f$'%(alpha_lin)) ),  loc='lower left', frameon=True,))
    
    ncol = 3
    axEvw_Sac.legend(loc=4, fontsize=FS+1, ncol=ncol, handlelength=2, columnspacing=1.2, labelspacing=0.3, fancybox=False, bbox_to_anchor=(1.05,-0.55))   
    
    #----------------------
    # Model config
    #----------------------
    
    props = dict(boxstyle='square', facecolor='wheat', alpha=0.5)
    textstr1 = '\n'.join(( r'"%s"'%(exprref.replace('_', '\_')), r'$L = %i$'%(L), r'$\nu = %.2e$'%(nu) ))
    axeigvals.text(-1.6, -0.1, textstr1, transform=axeigvals.transAxes, fontsize=FS, bbox=props)

    #----------------------
    # Save figure
    #----------------------
    fout = 'solutions/LATROT_%s__%i.png'%(exprref, tt+1)
    print('Saving %s'%(fout))
    plt.savefig(fout,dpi=dpi)

