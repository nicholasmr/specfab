# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
from matplotlib import ticker 
import matplotlib.cm as cm
import cmasher as cmr

sys.path.append('../')
from localheader import *

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import constants as sfconst
from specfabpy import statespace as sfsp
from specfabpy import plotting as sfplt

FS = sfplt.setfont_tex(fontsize=12)
FSANNO = FS
FSAX   = FS + 1
FSLBL  = FS + 2.5

norm = 1/np.sqrt(4*np.pi)

USE_LMAX_AS_REF = 1

#c_L10 = sfplt.c_dred
c_L10 = sfplt.c_dgreen

#--------------------
# Setup
#--------------------

L = 10 # paper assuems L=10, do not change
lm, nlm_len = sf.init(L) 

# Grain constants
Ecc, n_grain = 1, 1

# Calibration target of bulk behaviour for nlm_ref
Emt_ref = 8

#--------------------
# Integrate parcel for uniaxial compression
#--------------------

print('DYE 3 tests done on ice from z/H=%.2f to %.2f'%(1890/2037, 2006/2037))

Nt = 200 # time steps for Lagragian parcel integration

kwargs_LROT = dict(iota=1, Gamma0=None, nu=1) #      
DK = dict(type='ps', q=+0, tau=+1, axis=2)
strain_target = -0.95 # simulate parcel deformation until this target strain
nlm_uc, F_uc, time_uc, ugrad_uc = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, **kwargs_LROT)

DK = dict(type='ps', q=+0, tau=-1, axis=2)
strain_target = +6 # simulate parcel deformation until this target strain
nlm_ue, F_ue, time_ue, ugrad_ue = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, **kwargs_LROT)

#--------------------
# Determine maps
#--------------------

### Eij maps

nlm_unidir = sf__.nlm_ideal(m, 0, L__)
nlm_ref = nlm_uc[-1,:]
    
RES = 550 # production
#RES = 363 # faster
#RES = 100 # debug

Eca_list      = np.logspace(+0, np.log10(2000), int(RES)) # shear enhancement along basal plane
alpha_list    = np.linspace(0, 1, RES) # Sachs--Taylor weight 
logalpha_list = np.logspace(-2, 0, RES) # Sachs--Taylor weight 

#print(alpha_list)
#print(logalpha_list)

### Enhancement factors in state space

RESX = RESY = 8*50 
xlims = [-1.42,2.525]
ylims = [-1.5,3.55]

isvalid, x,y = sfsp.nlm_isvalid_grid(xlims, ylims, RESX, RESY)
imdat = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
for ii in range(RESY):
    for jj in range(RESX):
        imdat[ii,jj,:-1] = matplotlib.colors.to_rgb('w') if isvalid[ii,jj] else matplotlib.colors.to_rgb(sfplt.c_lgray)

#--------------------
# Plot
#--------------------

ii_Eca = np.argmin(np.abs(Eca_list - 10**3))
Eca = Eca_list[ii_Eca]
Emm_map, Emt_map, Epq_map, alpha_map, Eca_map = Eij_maps_tranisotropic(alpha_list, Eca_list, 'Mixed', Ecc=Ecc, n_grain=n_grain, nlm=nlm_ref)
alpha = alpha_map[0, np.argmin(np.abs(Emt_map[ii_Eca,:]-Emt_ref))]
Eij_grain = (Ecc,Eca)
grain_params = (Eij_grain, alpha, n_grain)

# Verify unidirectional limit reproduces calibration target for bulk behaviour (Emt)
Emt_actual = sf.Evw_tranisotropic(nlm_ref, m,t,tau_mt, *grain_params)
Emm_actual = sf.Evw_tranisotropic(nlm_ref, m,m,tau_mm, *grain_params)
print('Exz(Eca=%i, alpha=%.3f, nlm_ref) = %.3f (should be %.3f)'%(Eca, alpha, Emt_actual, Emt_ref))
print('Ezz(Eca=%i, alpha=%.3f, nlm_ref) = %.3f'%(Eca, alpha, Emm_actual))

### Setup figure

scale = 5.4
fig = plt.figure(figsize=(1.6*scale,0.7*scale))
gs = gridspec.GridSpec(1,2, width_ratios=[1.0,1.5])
gs.update(left=0.075, right=1-0.02, top=0.92, bottom=0.155, wspace=0.25, hspace=0.3)
ax1, ax2 = plt.subplot(gs[:, 1]), plt.subplot(gs[0, 0])

lw = 0.9
cmap = cm.RdBu
cmap.set_under('#fee6ce')
lvls = [0.8, 1, 2, 4, 6, 8, 10]
bnorm = matplotlib.colors.CenteredNorm(vcenter=1)

#--------------------
# Calibration map
#--------------------

def plot_calibmap(ax, nlm, alpha_list, pcalib,pmrk,pname):

    plt.sca(ax)
    logalpha = alpha_list[0] > 0
    if logalpha: ax.set_xscale('log') 
    ax.set_yscale('log')

    Emm_map, Emt_map, Epq_map, alpha_map, Eca_map = Eij_maps_tranisotropic(alpha_list, Eca_list, 'Mixed', Ecc=Ecc, n_grain=n_grain, nlm=nlm)

    hmap = ax.contourf(alpha_map, Eca_map, Emt_map, lvls, cmap=cmap, norm=bnorm, extend='both') # https://matplotlib.org/3.3.2/gallery/images_contours_and_fields/contourf_log.html

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
    
    if logalpha: logmid = (-0.9, 1)
    else:        logmid = (-0.1, 0.3)

    if logalpha: lvls_major = [0.5,]
    else:        lvls_major = [0.5,] #, 0.6, 0.7, 0.8, 0.9]
    kwargs = dict(linestyles='--',colors='k', linewidths=lw)
    CS2 = ax.contour(alpha_map, Eca_map, Emm_map, lvls_major, **kwargs)
    ax.clabel(CS2, CS2.levels, fmt='$E_{zz}=%.2f$', inline_spacing=10, manual=getlblpos(CS2, logmid))

    # plot second half here since contour labels need to be offset diffrently
    if logalpha: lvls_minor = [0.01, 0.02, 0.05, 0.1, 0.25]
    else:        lvls_minor = [0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
    CS2 = ax.contour(alpha_map, Eca_map, Emm_map, lvls_minor, **kwargs)
    ax.clabel(CS2, CS2.levels, fmt='$%.2f$', inline_spacing=10, manual=getlblpos(CS2, logmid))

    ### Misc

    hcalib, = ax.plot(pcalib[0], pcalib[1], pmrk, label=pname, markeredgecolor='k',markerfacecolor='w', markersize=11, clip_on=False, zorder=100)

    sfplt.panellabel(ax, 1, pname, frameon=True, pad=0.4, fontsize=FS)

    ax.set_xlabel(r'$\leftarrow$ Sachs \qquad\quad $\alpha$ \qquad\quad Taylor $\rightarrow$', fontsize=FSAX)
    ax.set_ylabel(r'$E_{ca}^\prime$', fontsize=FSAX)

    ax.set_xlim(alpha_list[[0,-1]])
    ax.set_ylim(Eca_list[[0,-1]])
    
    return hmap, hcalib


CS, hcalib = plot_calibmap(ax2, nlm_ref, alpha_list, (alpha, Eca), 'X', r'$\epsilon_{zz}=-0.95$')
hcb = plt.colorbar(CS, orientation='horizontal', pad=0.22, fraction=0.05, aspect=20)
hcb.set_label('$E_{xz}$, $E_{yz}$', fontsize=FSAX)
ax2.text(1.045, 1.6e3, r'Asymptotic $\rightarrow$', va='top', ha='center', rotation=90, fontsize=FSAX)
#ax2.legend([hcalib,], [r'Best fit: $E_{ca}^\prime,\, \alpha = %i,\, %.3f$'%(Eij_grain[1],alpha), ], bbox_to_anchor=(1,1.175), handletextpad=0.1, frameon=False)
ax2.text(alpha, Eca*0.7, 'Best fit', color='w', ha='center', va='top')
sfplt.panellabel(ax2, 2, r'\textbf{a}', frameon=False, bbox=(-0.24,1.175), fontsize=FSLBL)

#--------------------
# State map plot
#--------------------

ax = ax1
plt.sca(ax)

### Plot valid subspace

im = ax.imshow(imdat, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
ax.text(1.95, -1.1, '{\\bf Unphysical}\n\\bf{eigenvalues}', color=sfplt.c_dgray, ha='center', fontsize=FSANNO)

### Plot Eij

(Emm_smap, Emt_smap, x, y, isvalid) = sfsp.Eij_statemap_ice(grain_params, xlims, ylims, resx=RESX, resy=RESY)
CS = ax.contourf(x, y, Emt_smap, levels=lvls, extend='both', cmap=cmap, norm=bnorm, zorder=5)
#hcb = plt.colorbar(CS, orientation='horizontal', pad=0.17, fraction=0.05, aspect=26)
#hcb.set_label('$E_{xz}$, $E_{yz}$', fontsize=FSAX)

# determine locations for major (labelled) ticks 
lvls_major = np.array([0.4, ])
lvls_minor = np.array([0.2,0.6,0.8,1.0,1.3,1.6])
Ikeep = []
x0ref = 0.6
Ix = np.argmin(np.abs(x-x0ref))
for jj, lvl in enumerate(lvls_minor):
    Iy = np.nanargmin(np.abs(Emm_smap[:,Ix]-lvl))
    if -1.3 < y[Iy] < 1.9: Ikeep.append(jj)
manual_locations = [ (x0ref, y_) for y_ in np.linspace(-1.3,1.2,len(lvls_minor))[Ikeep] ]

kwargs = dict(linewidths=lw-0.2, linestyles='--', colors=sfplt.c_vdgray, zorder=10)

CS = ax.contour(x, y, Emm_smap, levels=lvls_minor, **kwargs)
ax.clabel(CS, CS.levels, inline=True, fmt=r'$%.1f$', fontsize=FS, manual=manual_locations)

CS  = ax.contour(x, y, Emm_smap, levels=lvls_major, **kwargs)
ax.clabel(CS, CS.levels, inline=True, fmt=r'$E_{zz}=%.1f$', fontsize=FS, manual=[(x0ref,1.8),])

### Misc

sfsp.plot_nlm_cases(ax, FSANNO, dx0_unidir=-0.06, dy0=0.06, isolbl_above=False, show_circle=False, c_planar=sfplt.c_vdgray, c_unidir=sfplt.c_vdgray, c_iso=sfplt.c_vdgray)

sfplt.panellabel(ax, 2, r'$E_{ca}^\prime,\, \alpha = %i,\, %.3f$'%(Eij_grain[1],alpha), frameon=True, pad=0.4, fontsize=FS)

plt.xlabel(r'${n}_2^0/n_0^0$', fontsize=FSAX)
plt.ylabel(r'${n}_4^0/n_0^0$', fontsize=FSAX)
plt.xlim(xlims)
plt.ylim(ylims)

#sfplt.panellabel(ax, 2, r'\textit{(c)}', frameon=False, bbox=(-0.13,1.1), fontsize=FSLBL)
sfplt.panellabel(ax, 2, r'\textbf{b}', frameon=False, bbox=(-0.13,1.12), fontsize=FSLBL)

### Parcel state trajectory

kwargs = dict(lw=1.75, c=c_L10, zorder=10) 
kwargs_mrk = dict(marker='o', ms=7.5, ls='none', zorder=20, c=c_L10)

n00, n20, n40 = nlm_uc[:,0], nlm_uc[:,sf.I20], nlm_uc[:,sf.I40]
nhat20, nhat40 = np.divide(n20,n00), np.divide(n40,n00)
ax.plot(nhat20, nhat40, '-', label='Uniaxial compression'%(), **kwargs)
ax.plot(nhat20[-1], nhat40[-1], **kwargs_mrk)

n00, n20, n40 = nlm_ue[:,0], nlm_ue[:,sf.I20], nlm_ue[:,sf.I40]
nhat20, nhat40 = np.divide(n20,n00), np.divide(n40,n00)
ax.plot(nhat20, nhat40, '--', label='Uniaxial extension'%(), **kwargs)
ax.plot(nhat20[-1], nhat40[-1], **kwargs_mrk)

legkwargs = {'frameon':False, 'ncol':1, 'fancybox':False, 'edgecolor':'k', 'labelspacing':0.3, 'handletextpad':0.7, 'handlelength':1.4}
hleg = ax.legend(loc=2, fontsize=FS, bbox_to_anchor=(0,0.88), **legkwargs)
hleg.get_frame().set_linewidth(0.75);


### Ideal CPO insets

if 1:

    geo, prj = sfplt.getprojection(rotation=45, inclination=50)

    W = H = 0.19 # ax width
    cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
    cmap.set_under('w')

    lvlmin, lvlmax = 0.2, 1.0
    ODF_plots = (\
        {'nlm':nlm_uc[-1,:], 'Fzz':F_uc[-1,2,2]-1, 'lvlmax':lvlmax, 'lvlmin':lvlmin}, \
        {'nlm':nlm_ue[-1,:], 'Fzz':F_ue[-1,2,2]-1, 'lvlmax':lvlmax, 'lvlmin':lvlmin}, \
    )

    for ODF in ODF_plots:
        lvlset = [np.linspace(ODF['lvlmin'], ODF['lvlmax'], 8), lambda x,p:'%.1f'%x]
        title = r'$\epsilon_{zz}=%.2f$'%(ODF['Fzz'])
        axin = sfplt.plotODFinset(ax,fig,prj, ODF['nlm'],lm, -90,0.8,0.13,'N', W,H, title=title, pad0mul=2.2, fstitle=FS, carr=c_L10, cmap=cmap, showcb=False, lvlset=lvlset)
        sfplt.plotcoordaxes(axin, geo, axislabels='vuxi', color=sfplt.c_dred,)
        

#--------------------
# Save figure
#--------------------

fname = 'Eij-calibrate.pdf'
print('Saving %s'%(fname))
plt.savefig(fname, transparent=0,  dpi=250)
plt.close()

