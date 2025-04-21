# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023-2024

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp
import pandas as pd
from progress.bar import Bar

OLIDATAPATH = '../../state-space/olivine'
sys.path.insert(0, OLIDATAPATH)
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import common as sfcom
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSSMALL = FS

import cartopy.crs as ccrs

### Flags

ENABLE_DDM = 0

#----------------------
# Experiment
#----------------------

T_EXP_UC = 1 # confined vertical compression 
T_EXP_SS = 2 # vertical shear

### Select experiment

for T_EXP in (T_EXP_UC, T_EXP_SS):

    #----------------------
    # Experiment definitions
    #----------------------

    Nt = 300
    L  = 10

    if T_EXP==T_EXP_SS:
        i,j = 0,2 # Fij components of interest
        DK = dict(type='ss', plane=1)
        strain_target = np.deg2rad(81) # simulate  parcel deformation until this target strain
        ODF_strains = [4, 4] if ENABLE_DDM else [1, 4] # show ODFs at these strains
        xlims = [0, 5.5]
        ylims = [1e-1, 4.5e0]
        expr = expr_DREX_ss # for DDM init state
                    
    if T_EXP==T_EXP_UC:
        i = j = 2 # Fij components of interest
        DK = dict(type='ps', q=0, axis=i)
        strain_target = -0.99 # simulate  parcel deformation until this target strain
        ODF_strains = [0.2, 0.2] if ENABLE_DDM else [0.2, 0.5] # show ODFs at these strains
        xlims = [-1, 0]
        ylims = [4e-2, 4.5e0]
        expr = expr_DREX_uc # for DDM init state
            
    def strainax(F):
        if T_EXP==T_EXP_SS: return [ 2*sf.F_to_strain(F[tt,:,:])[0,2] for tt in range(F.shape[0]) ] # https://www.continuummechanics.org/strain.html
        if T_EXP==T_EXP_UC: return [ sf.F_to_strain(F[ii,:,:])[2,2] for ii in range(F.shape[0]) ]
        
    #----------------------
    # CPO evolution
    #----------------------

    ### SDM

    lm, nlm_len = sf.init(L)
    blm = np.zeros((Nt+1,nlm_len), dtype=np.complex64) 
    nlm = np.zeros((Nt+1,nlm_len), dtype=np.complex64) 

    # Slip plane normal distribution
    kwargs_LROT = dict(iota=+1, Gamma0=None, nu=1)
    nlm[:,:], F, time, ugrad = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, **kwargs_LROT)

    # Slip direction distribution
    kwargs_LROT = dict(iota=-1, Gamma0=None, nu=1)
    blm[:,:], F, time, ugrad = sfint.lagrangianparcel(sf, DK, strain_target, Nt=Nt, **kwargs_LROT)

    ### DDM
                
    fname0 = expr['flist'][0] # initial state from DREX ensemble
    fname1 = '%s-%i.csv'%(fname0[:-4],1)
    fname2 = '%s-%i.csv'%(fname0[:-4],2)

    def load_axes(fullfname):
        df = pd.read_csv(fullfname, skiprows=expr['skiprows'], header=None, sep=expr['sep'])
        return df.to_numpy() # vj

    axes1_0 = load_axes('%s/data/%s/%s'%(OLIDATAPATH,expr['path'],fname1))
    axes2_0 = load_axes('%s/data/%s/%s'%(OLIDATAPATH,expr['path'],fname2))

    b0 = axes1_0 # grain no., xyz
    n0 = axes2_0
        
    N = Nt+1
    eps = np.array([(ugrad+ugrad.T)/2]*N)
    omg = np.array([(ugrad-ugrad.T)/2]*N)
    dt = dt = time[1]-time[0]
    
    if ENABLE_DDM:
        bi = sf.ri_LROT(b0, dt, N, eps,omg, -1) # time, grain no., xyz
        ni = sf.ri_LROT(n0, dt, N, eps,omg, +1)
        vi = np.array([np.cross(bi[nn,:,:], ni[nn,:,:]) for nn in range(N)])

    #----------------------
    # Determine eigenenhancements etc.
    #----------------------

    grain_params = sfconst.olivine['viscoplastic']['linear'] # Optimal n'=1 Sachs grain parameters
    print('Grain params: ', grain_params)

    ei, lami = sfcom.eigenframe(nlm, modelplane='xz')
    e1,e2,e3 = ei[:,0], ei[:,1], ei[:,2] # nn, i, xyz
    Eij   = sf.Eij_orthotropic_arr(nlm,blm,0*nlm, e1,e2,e3, *grain_params)

    x1,x2,x3 = sfcom.xi_tile(N)
    Exixj = sf.Eij_orthotropic_arr(nlm,blm,0*nlm, x1,x2,x3, *grain_params)
    
    if ENABLE_DDM:
        Eij_dsc = sf.Eij_orthotropic_discrete_arr(bi,ni,vi, e1,e2,e3, *grain_params)

    #----------------------
    # Plot
    #----------------------

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import matplotlib.gridspec as gridspec

    geo, prj = sfplt.getprojection(rotation=50+180, inclination=50)
    ODF_tsteps = np.array([np.argmin(np.abs(F[:,i,j]-thres)) for thres in ODF_strains])

    ### Setup figure

    scale = 0.87
    fig = plt.figure(figsize=(3.7*scale,3.7*scale), constrained_layout=False)
    gs = gridspec.GridSpec(1, 1)
    gs.update(hspace=0.53, wspace=0.26, left=0.185, right=0.95, top=0.92, bottom=0)
    ax_Y = fig.add_subplot(gs[0, 0])
    box = ax_Y.get_position()
    dy = box.height*0.5
    ax_Y.set_position([box.x0, box.y0+dy, box.width, box.height-dy])

    ### Plot relative strain-rates

    # Model relative lines
    
    lw = 1.5
    cshear = '#b35806'
    ccompr = '#542788'

    I11,I22,I33, I23,I13,I12 = 0,1,2, 3,4,5 # Voigt ordering
    Ixx,Iyy,Izz, Iyz,Ixz,Ixy = 0,1,2, 3,4,5 

    if ENABLE_DDM:
        # SI plot 
        ax_Y.semilogy(strainax(F), Eij[:,I12],     '-',  c='k', lw=lw, label=r"$E_{12}$", zorder=10)
        ax_Y.semilogy(strainax(F), Eij[:,I11],     '--', c='k', lw=lw, label=r"$E_{11}$", zorder=10)
        ax_Y.semilogy(strainax(F), Eij_dsc[:,I12], '-',  c=sfplt.c_green, label=r"$E_{12}$ DDM", zorder=9)
        ax_Y.semilogy(strainax(F), Eij_dsc[:,I11], '--', c=sfplt.c_green, label=r"$E_{11}$ DDM", zorder=9)
    else:
        # main text plot 
        ax_Y.semilogy(strainax(F), Exixj[:,Ixz], '-',  c=cshear, lw=lw, label=r"$E_{xz}$", zorder=10)
        ax_Y.semilogy(strainax(F), Exixj[:,Izz], '--', c=ccompr, lw=lw, label=r"$E_{zz}$", zorder=10)
        ax_Y.semilogy(strainax(F), Eij[:,I12], '-',  c='0.5', lw=lw, label=r"$E_{12}$", zorder=9)
        ax_Y.semilogy(strainax(F), Eij[:,I11], '--', c='0.5', lw=lw, label=r"$E_{11}$", zorder=9)

    # Isotropic line

    FSANNO = FS-1
    ciso = '0.0'
#    dyiso = 0.25
    if T_EXP==T_EXP_SS: xiso,dxiso,dyiso = 3.7+0.5, -1.6, 0.275
    if T_EXP==T_EXP_UC: xiso,dxiso,dyiso = 0.015-1, 0.275, 0.30
    ax_Y.plot(xlims, [1,1], ':', lw=1.2,  color=ciso, zorder=1)
    ax_Y.text(xiso+dxiso, 1, r'Isotropic', color=ciso, backgroundcolor='w', bbox=dict(boxstyle='square,pad=0.15', fc='w', ec='none'), ha='left', va='center', fontsize=FSANNO)
    ax_Y.text(xiso, 1+dyiso, r'$\uparrow$ Softer', c=ciso, ha='left', va='center', fontsize=FSANNO)
    ax_Y.text(xiso, 1-dyiso, r'$\downarrow$ Harder', c=ciso, ha='left', va='center', fontsize=FSANNO)

    # Axis labels

    F_indices = ('zz' if T_EXP==T_EXP_UC else 'xz')
    xaxlbl = r'\gamma_{%s}'%(F_indices) if T_EXP==T_EXP_SS else r'\epsilon_{%s}'%(F_indices)
    ax_Y.set_xlabel('$%s$'%(xaxlbl), labelpad=-0, fontsize=FS+1)
    ax_Y.set_ylabel(r'$E_{ij}$', labelpad=0)

    # Legend

    legkwargs = {'frameon':False, 'ncol':5, 'handlelength':1.05, 'handletextpad':0.25, 'columnspacing':0.6, 'fontsize':FSSMALL}
    leg = ax_Y.legend(loc=1, bbox_to_anchor=(1.05,1.25), **legkwargs)
    leg.get_frame().set_linewidth(0.8);

    # Show Hansen et al. obs limits?
    
    if T_EXP == T_EXP_SS and not ENABLE_DDM:
        strain = [4,10]
        Es = np.array([2e0, 4e0])
        Ec = np.array([2.5e-1, 9e-1])
        kwargsfill = {'ec':'none', 'alpha':1}
        ax_Y.fill_between(strain, [Es[0]]*2, [Es[1]]*2, color='#fee0b6', **kwargsfill)
        ax_Y.fill_between(strain, [Ec[0]]*2, [Ec[1]]*2, color='#d8daeb', **kwargsfill)
        ax_Y.annotate("Hansen et al. (2016)",
            xy=(4.35, 3.4e-1), xycoords='data',
            xytext=(2.2, 1.3e-1), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3,rad=-0.4", fc='k'))

    # Axis ticks and limits

    if T_EXP==T_EXP_SS: xticks = np.arange(0,20,0.5) 
    if T_EXP==T_EXP_UC: xticks = np.arange(-1,0.1,0.125)
    ax_Y.set_xticks(xticks[::2])  
    ax_Y.set_xticks(xticks[::], minor=True)  
    ax_Y.set_xlim(xlims)
    ax_Y.set_ylim(ylims)

    ### Plot ODFs

    for ii,nn in enumerate(ODF_tsteps):

        W = 0.17 # axis width
        fx = lambda ii, norb: 0.13 + ii*0.45 + norb*0.19 
        y = 0.01
        axb = plt.axes([fx(ii,0),y, W,W], projection=prj)
        axn = plt.axes([fx(ii,1),y, W,W], projection=prj)
 
        axb.set_global()
        axn.set_global()
        
        if ENABLE_DDM and ii==1:
            n_lat, n_colat, n_lon = sfdsc.cart2sph(ni[nn,:,:], deg=True)
            b_lat, b_colat, b_lon = sfdsc.cart2sph(bi[nn,:,:], deg=True)
            kwargs = dict(marker='o', s=0.6, linewidths=0.17, transform=geo, zorder=3)
            axn.scatter(n_lon, n_lat, c=sfplt.c_blue, edgecolors=sfplt.c_blue, facecolors=sfplt.c_blue, **kwargs)  
            axb.scatter(b_lon, b_lat, c=sfplt.c_red,  edgecolors=sfplt.c_red,  facecolors=sfplt.c_red,  **kwargs) 

            kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
            gl1 = axb.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
            gl2 = axn.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
            gl1.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
            gl2.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
        else:
            kwargs = dict(showcb=False, cbtickintvl=8, lvlset=[np.linspace(0,0.41,9), lambda x,p:'%.1f'%x], nchunk=None)        
            sfplt.plotODF(blm[nn,:], lm, axb, cblabel=r'$b/N$', cmap='Reds', **kwargs)
            sfplt.plotODF(nlm[nn,:], lm, axn, cblabel=r'$n/N$', cmap='Blues', **kwargs)

        axb.set_title(r'$b$', fontsize=FS)
        axn.set_title(r'$n$', fontsize=FS)
        sfplt.panellabel(axb, 2, r'$%s=%.1f$'%(xaxlbl, strainax(F)[nn]), bbox=(0.3,2.0), frameon=False, fontsize=FS)

        sfplt.plotcoordaxes(axb, geo, axislabels='vuxi', color='k')
        sfplt.plotcoordaxes(axn, geo, axislabels='vuxi', color='k')
        
        if not ENABLE_DDM:
        
            colors = ('#b3de69','#ffed6f','none')
            kwargs = dict(marker='.', ms=9, markeredgecolor=('k','k',None), markeredgewidth=0.5, colors=colors) # mec='k',
            sfplt.plotmi(axb, ei[nn], geo, **kwargs)
            sfplt.plotmi(axn, ei[nn], geo, **kwargs)

            if T_EXP==T_EXP_UC: x0,y0, dx,dy = (-1.17,4e-3, 0.03, 1.5e-3)
            if T_EXP==T_EXP_SS: x0,y0, dx,dy = (-0.9,1.6e-2, 0.03*5, 5e-3)
            kwargs = dict(marker='o', edgecolors='k', linewidths=0.5, s=18, clip_on=False, )
            ax_Y.scatter([x0,],[y0- 0,], c=colors[0], **kwargs)
            ax_Y.scatter([x0,],[y0-dy,], c=colors[1], **kwargs)
            ax_Y.text(x0+dx, y0- 0, r'$\vb{m}_1$', ha='left', va='center')
            ax_Y.text(x0+dx, y0-dy, r'$\vb{m}_2$', ha='left', va='center')

    ### Plot parcel deformation as inset

    if not ENABLE_DDM:

        if T_EXP==T_EXP_SS: parcel_tsteps = ODF_tsteps[[0,]]
        if T_EXP==T_EXP_UC: parcel_tsteps = ODF_tsteps[[1,]]

        for ii,nn in enumerate(parcel_tsteps):

            # Parcel locations (normalized figure coords, so must be set manually like this)
            if T_EXP==T_EXP_SS: 
                y0 = 0.46
                pc = np.array([[0.21,y0],])
                
            if T_EXP==T_EXP_UC: 
                y0 = 0.46
                pc = np.array([[0.5,y0]])

            axs = 0.15
            ax_sub = fig.add_axes([pc[ii,0], pc[ii,1], axs, axs], projection='3d')
            ax_sub.patch.set_alpha(0.0)
            lw = 0.7
            sfplt.plotparcel(ax_sub, F[nn,:,:], lw=lw, lwinit=lw, fonttex=True)
            
    ### Save figure

    T_EXP_STR = {T_EXP_UC:'UC', T_EXP_SS:'SS'}
    fname = 'Eij-lin-sachs--%s%s.pdf'%(T_EXP_STR[T_EXP], '-SI' if ENABLE_DDM else '')
    print('Saving output to %s'%(fname))
    fig.patch.set_alpha(0.0)
    plt.savefig(fname, dpi=175)

