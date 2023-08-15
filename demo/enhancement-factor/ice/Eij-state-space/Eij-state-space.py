# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
from matplotlib import ticker 
import matplotlib.cm as cm
import cmasher as cmr

sys.path.append('../../../')
import demolib as dl

sys.path.append('../..')
from localheader import *

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt

FS = sfplt.setfont_tex(fontsize=12)
FSANNO = FS - 1
FSLEG  = FSANNO + 0.5
FSAX   = FSANNO + 2

norm = 1/np.sqrt(4*np.pi)

#--------------------
# Setup
#--------------------

# Grain constants
n_grain = 1
Ecc = 1

# Calibration target of bulk behaviour
Emt_ref = 10 

#--------------------
# Determine maps
#--------------------

### Eij for unidir

L = 4
lm, nlm_len = sf.init(L) 
nlm_unidir = sf__.nlm_ideal(m, 0, L__)
    
RES = 20 * 10
Eca_list   = np.logspace(+0, 3, int(RES*2)) # shear enhancement along basal plane
alpha_list = np.logspace(-2, 0, RES) # Sachs--Taylor weight 
Emm_map, Emt_map, Epq_map, alpha_map, Eca_map = Eij_maps_tranisotropic(alpha_list, Eca_list, 'Mixed', Ecc=Ecc, n_grain=n_grain)
Zpq = np.ma.array(np.divide(Emt_map,Epq_map))        

### Enhancement factors in state space

RESX = RESY = 8*50 
xlims = [-1.42,2.525]
ylims = [-1.5,3.55]

isvalid, x,y = dl.nlm_isvalid_grid(xlims, ylims, RESX, RESY)
imdat = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
for ii in range(RESY):
    for jj in range(RESX):
        imdat[ii,jj,:-1] = matplotlib.colors.to_rgb('w') if isvalid[ii,jj] else matplotlib.colors.to_rgb(sfplt.c_lgray)

#--------------------
# Plot
#--------------------

alpha0 = alpha_map[0, np.argmin(np.abs(Emt_map[-1,:]-Emt_ref))]
alpha_sweep = np.logspace(np.log10(alpha0), 0, 20)

for gg, alpha in enumerate(alpha_sweep):

    ii = np.argmin(np.abs(alpha-alpha_list))
    Eca = Eca_list[ np.argmin(np.abs(Emt_map[:,ii]-Emt_ref)) ]
    Eij_grain = (Ecc,Eca)
    grain_params = (Eij_grain, alpha, n_grain)
    
    # Verify unidirectional limit reproduces calibration target for bulk behaviour (Emt)
    Emt_actual = sf.Evw_tranisotropic(nlm_unidir, m,t,tau_mt, *grain_params)
    print('Emt(Eca=%i, alpha=%.3f) = %.2f (should be %.2f)'%(Eca, alpha, Emt_actual, Emt_ref))
    
    ### Setup figure

    scale = 5.1
    fig = plt.figure(figsize=(1.4*scale,0.8*scale))
    gs = gridspec.GridSpec(1,2, width_ratios=[1,1.5])
    gs.update(left=0.09, right=1-0.02, top=0.87, bottom=0.03, wspace=0.24)
    ax1, ax2 = plt.subplot(gs[0, 0]), plt.subplot(gs[0, 1])

    lw = 0.9

    ax1.set_xscale('log') 
    ax1.set_yscale('log')

    X, Y = alpha_map, Eca_map

    #--------------------
    # Plot 1
    #--------------------

    ax = ax1
    plt.sca(ax)
    
    cmap = cm.Blues
    lvls = [1, 2, 5, 10, 20]
    bnorm = matplotlib.colors.BoundaryNorm(lvls, 256, extend='both')
    hmap = ax.contourf(X, Y, Emt_map, lvls, cmap=cmap, norm=bnorm, extend='both') # https://matplotlib.org/3.3.2/gallery/images_contours_and_fields/contourf_log.html

    lvls_major = [0.02,  0.05, 0.1, 0.25, 0.5]
    CS2 = ax.contour(X, Y, Emm_map, lvls_major, linestyles='--',colors='k', linewidths=lw)

    ### Contour labels and their positions

    def getlblpos(CS,logmid):
        label_pos = []
        for line in CS.collections:
            for path in line.get_paths():
                logvert = np.log10(path.vertices)
                # find closest point
                logdist = np.linalg.norm(logvert-logmid, ord=2, axis=1)
                min_ind = np.argmin(logdist)
                label_pos.append(10**logvert[min_ind,:])
        return label_pos
        
    xmin,xmax,ymin,ymax = plt.axis()
    logmid_2 = (-1.1, 1)
    ax.clabel(CS2, CS2.levels, fmt='$E_{mm}=%.2f$', inline_spacing=10, manual=getlblpos(CS2, logmid_2))

    ### Misc

    hsim, = ax.plot(alpha, Eca, 'X', markeredgecolor='k',markerfacecolor='w', markersize=12, clip_on=False, zorder=100)

    hcb = plt.colorbar(hmap, orientation='horizontal', pad=0.18, aspect=18)
    hcb.set_label(r'$E_{mt}$', fontsize=FSAX)
    ax.set_xlabel(r'$\alpha$', fontsize=FSAX)
    ax.set_ylabel(r'$E_{ca}^\prime$', fontsize=FSAX)

    ax.set_xlim(alpha_list[[0,-1]])
    ax.set_ylim(Eca_list[[0,-1]])

    sfplt.panellabel(ax, 1, r'$(E_{cc}^\prime,n^\prime)=(%i,%i)$'%(Ecc,n_grain), fontsize=FS)
    ax.set_title('Calibration of grain parameters\n for unidirectional CPO', fontsize=FSAX, pad=10)

    #--------------------
    # Plot 2
    #--------------------

    ax = ax2
    plt.sca(ax)

    ### Plot valid subspace

    im = ax.imshow(imdat, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
    ax.text(1.95, -1.1, '{\\bf Unphysical}\n\\bf{eigenvalues}', color=sfplt.c_dgray, ha='center', fontsize=FSANNO)

    ### Plot Eij

    (Emm_smap, Emt_smap, x, y, isvalid) = dl.Eij_statemap_ice(grain_params, xlims, ylims, resx=RESX, resy=RESY)

    cmap = cm.RdBu
    cmap.set_under('#fee6ce')
    CS = ax.contourf(x, y, Emt_smap, levels=[0.5, 1, 1.5, 2,3,4,5], extend='both', cmap=cmap, norm=matplotlib.colors.CenteredNorm(vcenter=1), zorder=5)
    hcb = plt.colorbar(CS, orientation='horizontal', pad=0.18, aspect=25.5)
    hcb.set_label('$E_{mt}$', fontsize=FSAX)

    # determine locations for major (labelled) ticks 
    lvls_major = np.array([0.5, 1, 1.5])
    Ikeep = []
    x0ref = 0.6
    Ix = np.argmin(np.abs(x-x0ref))
    for jj, lvl in enumerate(lvls_major):
        Iy = np.nanargmin(np.abs(Emm_smap[:,Ix]-lvl))
        if -1.0 < y[Iy] < 1.9: Ikeep.append(jj)
    lvls_major = lvls_major[Ikeep]
    manual_locations = [(x0ref, y_) for y_ in np.array([1.2,0.2,-1.0])[Ikeep]]

    kwargs = dict(linewidths=lw-0.2, linestyles='--', colors=sfplt.c_vdgray, zorder=10)
    CS_ = ax.contour(x, y, Emm_smap, levels=np.array([0.25,0.75,1.25,1.75]), **kwargs)
    CS  = ax.contour(x, y, Emm_smap, levels=lvls_major, **kwargs)
    ax.clabel(CS, CS.levels, inline=True, fmt=r'$E_{mm}=%.1f$', fontsize=FSAX, manual=manual_locations)

    ### Misc

    dl.plot_nlm_cases(ax, FSANNO, dx0_unidir=-0.09, isolbl_above=False, show_circle=False)

    sfplt.panellabel(ax, 2, r'$(\alpha, E_{ca}^\prime) = (%.3f, %i)$'%(alpha, Eij_grain[1]), fontsize=FS)

    legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'columnspacing': 0.5, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}
    leg = plt.legend(loc=2, fontsize=FSLEG, frameon=False, ncol=1, **legkwargs); 
    plt.xlabel(r'$\hat{n}_2^0$', fontsize=FSAX)
    plt.ylabel(r'$\hat{n}_4^0$', fontsize=FSAX)
    plt.xlim(xlims)
    plt.ylim(ylims)

    ax.set_title('Eigenenhancements in CPO state space', fontsize=FSAX, pad=10)

    ### Ideal CPO insets
    
    if 1:

        geo, prj = sfplt.getprojection(rotation=55+180, inclination=50)

        arr = lambda ang: 0.125*np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
        n00 = norm
        lvlmax = 1.0
        lvlmin = 0.2
        ODF_plots = (\
            {'nlm':sf.nlm_ideal([0,0,1], 0, 8),       'title':'', 'axloc':(0.877, 0.585), 'darr':arr(-90), 'lvlmax':lvlmax, 'lvlmin':lvlmin*1.8}, \
            {'nlm':sf.nlm_ideal([0,0,1], np.pi/2, 8), 'title':'', 'axloc':(0.472, 0.375), 'darr':arr(-90), 'lvlmax':lvlmax, 'lvlmin':lvlmin}, \
        )

        for ODF in ODF_plots:

            W = 0.135 # ax width
            axpos = [ODF['axloc'][0],ODF['axloc'][1], W,W]
            axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
            axin.set_global()
            
            nlm = ODF['nlm']
            cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
            cmap.set_under('w')
            lvlset = [np.linspace(ODF['lvlmin'], ODF['lvlmax'], 8), lambda x,p:'%.1f'%x]
            sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False)

            # Arrow to ODF state
            n20_, n40_ = np.real(nlm[3])/norm, np.real(nlm[10])/norm
            sc = np.diff(ylims)/np.diff(xlims)
            ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                            xytext=(n20_+ODF['darr'][0]/norm, n40_+sc**2*ODF['darr'][1]/norm), textcoords='data', \
                            arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", linewidth=1.5, edgecolor='0.2', facecolor='0.2'),zorder=20)            

            axin.set_title(ODF['title'], fontsize=FS)

    #--------------------
    # Save figure
    #--------------------

    fname = 'frames/Eij-state-space-%i.png'%(gg)
    print('Saving %s'%(fname))
    plt.savefig(fname, transparent=0,  dpi=250)
    plt.close()

