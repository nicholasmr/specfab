# N. M. Rathmann <rathmann@nbi.ku.dk>, 2020-2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp
from netCDF4 import Dataset

from specfabpy import specfab as sf
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

expname = sys.argv[1]

#------------------
# Load solution
#------------------

fh = Dataset('solutions/LROT_%s.nc'%(expname), mode='r')
loadvar = lambda field: np.array(fh.variables[field][:])

# Model config
Nt, dt, L = fh.getncattr('tsteps'), fh.getncattr('dt'), fh.getncattr('L')

# Grain parameters
Eca_lin,  Ecc_lin  = fh.getncattr('Eca_lin'),  fh.getncattr('Ecc_lin')
Eca_nlin, Ecc_nlin = fh.getncattr('Eca_nlin'), fh.getncattr('Ecc_nlin')
alpha_lin, alpha_nlin = fh.getncattr('alpha_lin'), fh.getncattr('alpha_nlin')

# CPO state
lm, c = loadvar('lm'), loadvar('c_re') + 1j*loadvar('c_im') 
lm = np.array(lm).T
eigvals = loadvar('eigvals')
m1,m2,m3 = loadvar('m1'), loadvar('m2'), loadvar('m3')
p1,p2,p3 = loadvar('p1'), loadvar('p2'), loadvar('p3')

# Enhancement factors
Eij_lin, Eij_nlin = loadvar('Eij_lin'), loadvar('Eij_nlin') 
Epij_lin, Epij_nlin = loadvar('Epij_lin'), loadvar('Epij_nlin') 

#------------------
# Plot
#------------------

tsteps = [0, int(Nt*1/2), Nt-1] # time steps to plot
#tsteps = [Nt-1,] # time steps to plot

lw0, lw1, lw2 = 2.5,2.25,2.0

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from matplotlib.offsetbox import AnchoredText
    
rotation=55 
inclination=45
geo, prj = sfplt.getprojection(rotation=rotation, inclination=inclination)

def plot_vec(ax, v, lbl, color, ls='-', lw=2):
    ax.plot([0, +v[0]],[0, +v[1]],[0,+v[2]], color=color, ls=ls, lw=lw, label=lbl)
    ax.plot([0, -v[0]],[0, -v[1]],[0,-v[2]], color=color, ls=ls, lw=lw)

for tt in tsteps:

    #----------------------
    # Figure setup
    #----------------------

    dpi, scale = 200, 3.3
    fig = plt.figure(figsize=(3/2*2.1*scale,2.3*scale))
    gs = gridspec.GridSpec(2,3, height_ratios=[1,1.2], width_ratios=[1,1,1])
    gs.update(left=-0.03, right=1-0.06/3, top=0.97, bottom=0.20, wspace=0.015*18, hspace=0.35)

    ax_ODF     = plt.subplot(gs[0, 0], projection=prj)
    ax_mi      = plt.subplot(gs[0, 1], projection='3d')
    ax_eigvals = plt.subplot(gs[0, 2])
    ax_Elin    = plt.subplot(gs[1, 1])
    ax_Enlin   = plt.subplot(gs[1, 2])

    ax_ODF.set_global() # be sure to show entire S^2
                
    #----------------------
    # ODF
    #----------------------

    ax = ax_ODF
    sfplt.plotODF(c[tt,:], lm, ax, cmap='Greys')
    sfplt.plotcoordaxes(ax, geo, axislabels='vuxi')

    #----------------------
    # Eigenvalues
    #----------------------

    ax_eigvals.plot([tt,tt],[0,1],':k', lw=2)
    
    steps = np.arange(len(eigvals[:,0]))
    ax_eigvals.plot(steps,eigvals[:,0], '-', c=sfplt.c_red,   label='$a_{1}$', lw=lw0)
    ax_eigvals.plot(steps,eigvals[:,1], '-', c=sfplt.c_green, label='$a_{2}$', lw=lw1)
    ax_eigvals.plot(steps,eigvals[:,2], '-', c=sfplt.c_blue,  label='$a_{3}$', lw=lw2)
    
    ax_eigvals.set_ylim([0,1])
    ax_eigvals.set_xlim([0, Nt+1])
    ax_eigvals.set_xlabel('time step')
    ax_eigvals.set_ylabel('$a_{i}$')
    ax_eigvals.grid()              
    ax_eigvals.legend(handlelength=1, ncol=1, labelspacing=0.3, fancybox=False, loc=2)
            
    #----------------------
    # Principal frame
    #----------------------

    ax_mi.view_init(elev=90-inclination, azim=rotation) # same as ODF plot
    
    ax_mi.set_xlabel('$x$'); ax_mi.set_xlim([-1,1])
    ax_mi.set_ylabel('$y$'); ax_mi.set_ylim([-1,1])
    ax_mi.set_zlabel('$z$'); ax_mi.set_zlim([0,1])
    
    plot_vec(ax_mi,m1[tt,:], r'$\vb{m}_1$', sfplt.c_red)
    plot_vec(ax_mi,m2[tt,:], r'$\vb{m}_2$', sfplt.c_green)
    plot_vec(ax_mi,m3[tt,:], r'$\vb{m}_3$', sfplt.c_blue)
    
    lwpq=1
    plot_vec(ax_mi,p1[tt,:], r'$\vb{p}_{1}$', sfplt.c_red,   ls='--', lw=lwpq)
    plot_vec(ax_mi,p2[tt,:], r'$\vb{p}_{2}$', sfplt.c_green, ls='--', lw=lwpq)
    plot_vec(ax_mi,p3[tt,:], r'$\vb{p}_{3}$', sfplt.c_blue,  ls='--', lw=lwpq)
    
    ax_mi.legend(handlelength=1, loc=1, bbox_to_anchor=(1.17,1), fancybox=False)

    #----------------------
    # Enhancement factors
    #----------------------
    
    lblm = lambda ii,jj: '$E_{m_%i m_%i}$'%(ii+1,jj+1)
    lblp = lambda ii,jj: '$E_{p_%i p_%i}$'%(ii+1,jj+1)

    def plot_enhancements(ax, Eij, Epij):
        
        ax.semilogy(steps, Eij[:,0], '-', c=sfplt.c_red,   label=lblm(0,0), lw=lw0)
        ax.semilogy(steps, Eij[:,1], '-', c=sfplt.c_green, label=lblm(1,1), lw=lw1)
        ax.semilogy(steps, Eij[:,2], '-', c=sfplt.c_blue,  label=lblm(2,2), lw=lw2)    
        
        ax.semilogy(steps, Eij[:,3], ':', c=sfplt.c_lgreen, label=lblm(1,2), lw=lw2)
        ax.semilogy(steps, Eij[:,4], ':', c=sfplt.c_lblue,  label=lblm(0,2), lw=lw1)
        ax.semilogy(steps, Eij[:,5], ':', c=sfplt.c_lred,   label=lblm(0,1), lw=lw0)

        Eratio = np.divide(Eij[:,5], Epij[:,5])
        ax.semilogy(steps, Epij[:,5], '-', c=sfplt.c_gray, label=lblp(0,1), lw=lw2) 
        ax.semilogy(steps, Eratio, '--', c=sfplt.c_gray, label=lblm(0,1)+'/'+lblp(0,1), lw=lw2)

        xlims=[0, Nt+1]
        ax.semilogy(xlims, [2.5,2.5], '--k', lw=1) 
        ax.semilogy(xlims, [4.375,4.375], '--k', lw=1)
        ax.set_xlim(xlims)
        ax.set_ylim([np.amin([1e-1, np.amin(Eij[:]), np.amin(Epij[:])]), np.amax([1e+1, np.amax(Eij[:]), np.amax(Epij[:]), np.amax(Eratio)])])
        ax.set_xlabel('time step')
        ax.set_ylabel('$E_{vw}$')
        ax.grid()
    
    plot_enhancements(ax_Enlin, Eij_nlin, Epij_nlin)
    ax_Enlin.set_title(r'Nonlinear Sachs')
    
    plot_enhancements(ax_Elin,  Eij_lin,  Epij_lin)
    ax_Elin.set_title(r'Linear mixed Taylor--Sachs')
    
    ax_Enlin.add_artist(AnchoredText('\n'.join( (r"$n' = %i$"%(3), r'$E_{cc} = %.1f$'%(Ecc_nlin), r'$E_{ca} = %.0e$'%(Eca_nlin), r'$\alpha = %.4f$'%(alpha_nlin)) ), loc='lower left', frameon=True,))
    ax_Elin.add_artist( AnchoredText('\n'.join( (r"$n' = %i$"%(1), r'$E_{cc} = %.1f$'%(Ecc_lin),  r'$E_{ca} = %.0e$'%(Eca_lin),  r'$\alpha = %.4f$'%(alpha_lin)) ),  loc='lower left', frameon=True,))
    
    ax_Enlin.legend(loc=4, fontsize=FS+1, ncol=3, handlelength=2, columnspacing=1.2, labelspacing=0.3, fancybox=False, bbox_to_anchor=(1.05,-0.55))   
    
    #----------------------
    # Model config
    #----------------------
    
    props = dict(boxstyle='square', facecolor='wheat', alpha=0.5)
    textstr1 = '\n'.join(( r'"%s"'%(expname.replace('_', '\_')), r'$L = %i$'%(L) ))
    ax_eigvals.text(-1.6, -0.1, textstr1, transform=ax_eigvals.transAxes, fontsize=FS, bbox=props)

    #----------------------
    # Save figure
    #----------------------
    
    fout = 'solutions/LROT_%s__%i.png'%(expname, tt+1)
    print('Saving %s'%(fout))
    plt.savefig(fout, dpi=dpi)

