# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2024

"""
CPO state-space diagram comparing model trajectories to data
"""

# works with pip3 install matplotlib==3.7

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import pickle
import scipy

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colormaps as cm
from matplotlib.collections import LineCollection
import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["k", "r", "c"]) # black first for making legend of linecollections black
import cmcrameri.cm as cmc
import cartopy.crs as ccrs

sys.path.insert(0, '..')
from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import integrator as sfint
from specfabpy import statespace as sfsp
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-0.75
FSANNO = FS-1

SELFNAME = sys.argv[0][:-3] # used as prefix for pickled files
os.system('mkdir -p specfab-state-trajectories')
def pfile(fname): return "specfab-state-trajectories/%s--%s.p"%(SELFNAME, fname) # full path name for pickled files

#--------------------
# Run options
#--------------------

DEBUG           = 0  # Low-resolution plotting, etc.

INTEGRATE_MODEL = 0  # Generate model lines? Else load saved Pickle files.
PLOT_ODF_INSETS = 1
            
#--------------------
# Config
#--------------------

Nt = 200 # Number of specfab integration steps

RESX = RESY = 200 if DEBUG else 500

L = 20
lm, nlm_len = sf.init(L) 
    
#--------------------
# Experimental data to plot
#--------------------

#experiments_b = (expr_Yabe2020_b, expr_Miyazaki2013_uc_b, expr_Miyazaki2013_ue_b, expr_Kim2022_b, expr_Bernard2019_b)
#experiments_n = (expr_Yabe2020_n, expr_Miyazaki2013_uc_n, expr_Miyazaki2013_ue_n, expr_Kim2022_n, expr_Bernard2019_n)

#experiments = (expr_Boneh2014, expr_Kumamoto2019, expr_Bernard2019)
experiments = (expr_Boneh2014, expr_Yabe2020, expr_Miyazaki2013_uc, expr_Miyazaki2013_ue, expr_Kim2022) 
#experiments = (expr_Boneh2014, expr_Yabe2020)

FTYPE = 'A' # fabric type

mrk = ['d','s','o','o', 'X'] # should be same as in J-index-validation.py plot

#--------------------
# Modelled correlations
#--------------------

def integrate_model(nlm0, modtype, Mtype='LROT', Nt=Nt):
    iota = zeta = Gamma0 = Lambda = nu = None
    if Mtype == 'LROT': zeta, nu = 0, 1
    if modtype == 'uc': mod, target = dict(type='ps', axis=2, T=+1, r=0), -0.99
    if modtype == 'ue': mod, target = dict(type='ps', axis=2, T=-1, r=0), 10
    kwargs = dict(Nt=Nt, nlm0=nlm0, zeta=zeta, Gamma0=Gamma0, Lambda=Lambda, nu=nu)
    nlm, F, *_ = sfint.lagrangianparcel(sf, mod, target, iota=+1, **kwargs)
    blm, F, *_ = sfint.lagrangianparcel(sf, mod, target, iota=-1, **kwargs)
    strain_ij = np.array([sf.F_to_strain(F[nn,:,:]) for nn in np.arange(Nt)]) # strain_ij tensor
    strainzz = strain_ij[:,2,2] # vertical strain
    pickle.dump([blm, nlm, strainzz, lm, nlm_len], open(pfile(modtype), "wb"))
    #return blm, nlm, strainzz
    
    
if INTEGRATE_MODEL:
    # Solve for state-vector time evolution
    print('*** Generating model trajectories from scratch. You can re-use the saves trajectories to avoid re-calculating them by setting INTEGRATE_MODEL=1')
    nlm_iso = np.zeros((nlm_len), dtype=np.complex64)
    nlm_iso[0] = norm
    integrate_model(nlm_iso, 'uc') 
    integrate_model(nlm_iso, 'ue')
    
### Load solutions

def load_solution(expname):
    blm, nlm, strainzz, lm, nlm_len = pickle.load(open(pfile(expname), "rb"))    
    nlm  = np.array([ nlm[tt,:]/nlm[tt,0] for tt in np.arange(Nt+1) ]) # normalize
    blm  = np.array([ blm[tt,:]/blm[tt,0] for tt in np.arange(Nt+1) ]) # normalize
    return (blm, nlm, strainzz, lm, nlm_len)

blm_uc, nlm_uc, strainzz_uc, lm, nlm_len = load_solution('uc')
blm_ue, nlm_ue, strainzz_ue, lm, nlm_len = load_solution('ue')

#--------------------
# Plot
#--------------------

mse = 7.5 # end-member case points
ms  = mse # data points
mew = 1.5 # marker edge width

def plot_trajectory(ax, xlm, strainzz, cmap, cnorm, lw=2, ls='-'):
    x, y = np.real(xlm[:,sf.I20]), np.real(xlm[:,sf.I40])
    if ls != '-':
        f = scipy.interpolate.interp1d(x,y, bounds_error=False)
        xnew = np.linspace(0,-1.2, 20) if np.amin(x) < -0.1 else np.linspace(0, 2.3, 35)
        ynew = f(xnew)
        fs = scipy.interpolate.interp1d(x,strainzz, bounds_error=False)
        strainzz_new = fs(xnew)
        x, y, strainzz = xnew, ynew, strainzz_new

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=cnorm, ls=ls)
    lc.set_array(strainzz)
    lc.set_linewidth(lw)
    line = ax.add_collection(lc)
    return line

rng = [0,]
rng = range(4)
for ii in rng:

    scale = 4
    fig = plt.figure(figsize=(0.95*scale,0.7*scale))
#    ax = plt.gca()
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.13, bottom=0.19, right=0.99, top=0.98, wspace=0, hspace=0)
    
#    ax.patch.set_facecolor('white')
    pp1 = plt.Rectangle((-2,-2), 10, 10, color='w') 
    ax.add_patch(pp1) 
    
    axname = 'b' if ii in (0,2) else 'n'
    IS_UC = ii in (0,1)
    IS_UE = ii in (2,3)
    cmap = cmc.batlow_r if IS_UC else cmc.glasgow 
    cmax_ue = 3
    cnorm = plt.Normalize(-1, 0) if IS_UC else plt.Normalize(0, cmax_ue)
    
    xlims, ylims = [-1.35,2.55], [-1.35,3.5]
    sc = np.diff(ylims)/np.diff(xlims)

    ### Determine valid subspace

    isvalid, x,y = sfsp.nlm_isvalid_grid(xlims, ylims, RESX, RESY)
    C = statespace_shading(nlm_uc, nlm_ue, x, y, isvalid, shade_circle=False)
    im = ax.imshow(C, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)

    #--------------------
    # Construct plot
    #--------------------

    ### Model lines

    hwmul=0.75
    lw = 2
#    n20arrpos = 0.175
    
    # UC experiment
    if ii==0: xlm_sdm, strainzz_sdm = blm_uc, strainzz_uc
    if ii==1: xlm_sdm, strainzz_sdm = nlm_uc, strainzz_uc

    # UE experiment
    if ii==2: xlm_sdm, strainzz_sdm = blm_ue, strainzz_ue
    if ii==3: xlm_sdm, strainzz_sdm = nlm_ue, strainzz_ue

###

    line_sdm = plot_trajectory(ax, xlm_sdm, strainzz_sdm, cmap, cnorm)
    
    cax = fig.add_axes([0.845 + (0 if IS_UC else 0.02), 0.225, 0.020, 0.32])
    plt.sca(ax)
    ticks = np.arange(-1, 0+1e-3, 0.2) if IS_UC else np.arange(0, cmax_ue+1e-3, 0.5)
    cbar = plt.colorbar(line_sdm, cax=cax,  extend='neither' if IS_UC else 'max', ticks=ticks)
    cbar.ax.set_title('$\epsilon_{zz}$', fontsize=FS+2)

    ### D-Rex model lines
    
    drexexp = 'axisymmetricCompression' if IS_UC else 'axisymmetricExtension'
    ri = 1 if axname=='b' else 2 # crystallographic axis number
    xlm_drex, F = load_drex_run(drexexp, ri, normalize=True)

    strain_ij = np.array([sf.F_to_strain(F[nn,:,:]) for nn in np.arange(len(xlm_drex))]) # strain_ij tensor
    strainzz_drex = strain_ij[:,2,2]
    lines_drex = plot_trajectory(ax, xlm_drex, strainzz_drex, cmap, cnorm, lw=3, ls='dotted')

    ### Labels
        
    sfsp.plot_nlm_cases(ax, FSANNO, ms=7.5, show_circle=False, dy0=0.075, dx0_unidir=-0.09, dx0_planar=0.02, isolbl_above=True)

    kwargs_lbl = dict(ha='center', fontsize=FSANNO)
    ax.text(+0.300/norm, +0.150/norm, '{\\bf Single maximum}', color=sfsp.c_unidir, rotation=40, **kwargs_lbl)
    ax.text(-0.170/norm, +0.11/norm, '{\\bf Girdle}', color=sfsp.c_planar, rotation=-40, **kwargs_lbl)
#    ax.text(-0.250/norm, -0.320/norm, '{\\bf Unphysical}\n{\\bf eigenvalues}', color=sfplt.c_dgray, rotation=0, **kwargs_lbl)
    
    legkwargs = {'handlelength':1.6, 'framealpha':1.0, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}
    leg_model = ax.legend([line_sdm,lines_drex], ['SDM','D-Rex'], loc=2, bbox_to_anchor=(0.0,1.03), ncol=2, columnspacing=1.3, fontsize=FSLEG, frameon=False, **legkwargs)
    
    ### Experimental data points

    correlations = []

    for expr in experiments:
        fcorr = "observed-states/%s.p"%(expr['path']) 
        print('=== Loading correlations: %s ==='%(fcorr))
        corr = pickle.load(open(fcorr, "rb"))
        corr['%s20'%(axname)] /= norm # Normalization
        corr['%s40'%(axname)] /= norm
        corr['%s60'%(axname)] /= norm
        correlations.append(corr) 

    handles = []
    titles  = []

    for jj,expr in enumerate(experiments):
        if IS_UC and expr['type'] != 'uc': continue
        if IS_UE and expr['type'] != 'ue': continue
        
        x = correlations[jj]['%s20'%(axname)]
        y = correlations[jj]['%s40'%(axname)]
#        if expr['type']=='ss':  mrk = 's'
#        if expr['type']=='ue':  mrk = '^'
#        if expr['type']=='uc':  mrk = 'd'
#        if expr['type']=='ucw': mrk = 's'
#        if expr['type']=='cc':  mrk = 'X'
        fillstyle = 'none' if expr['path'] != 'Boneh_etal_2014' else 'full'
        strainzz = np.array(expr['Fzz']) - 1
        for kk in range(len(x)):
            c = cmap( (strainzz[kk]-np.amin(ticks))/(np.amax(ticks)-np.amin(ticks)) )
            h, = ax.plot(x[kk],y[kk], ms=9, markeredgewidth=mew, ls='none', color=c, fillstyle=fillstyle, marker=mrk[jj])
        handles.append(h)
        titles.append(expr['plotname'])

    leg_obs = ax.legend(handles, titles, loc=2, bbox_to_anchor=(0.0,0.94), ncol=1, fontsize=FSLEG, frameon=False, **legkwargs)
    [lgd.set_color('black') for lgd in leg_obs.legendHandles]
    ax.add_artist(leg_model)
    
    ### Plot ODF insets?

    if PLOT_ODF_INSETS:

        ### Define model states
        
        arr = lambda ang: 0.1*np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
        strainzz_target = -0.5 if IS_UC else 2
        ttl_sdm  = 'SDM\n $\epsilon_{zz}=%.1f$'%(strainzz_target)
        ttl_drex = 'D-Rex\n $\epsilon_{zz}=%.1f$'%(strainzz_target)
        axloc_sdm  = (0.46, 0.11)
        axloc_drex = (axloc_sdm[0]+0.20, axloc_sdm[1])
        xlm_sdm_  = xlm_sdm[ np.argmin(np.abs(strainzz_sdm -strainzz_target)),:]
        xlm_drex_ = xlm_drex[np.argmin(np.abs(strainzz_drex-strainzz_target)),:]
        cmap = 'Reds' if axname=='b' else 'Blues'
        lvlmax = 0.3 if IS_UC else 0.4
        ODF_plots = (\
            {'xlm':xlm_sdm_,  'title':ttl_sdm,  'cax':None, 'cmap':cmap, 'axloc':axloc_sdm,  'lvlmax':lvlmax}, \
            {'xlm':xlm_drex_, 'title':ttl_drex, 'cax':None, 'cmap':cmap, 'axloc':axloc_drex, 'lvlmax':lvlmax}, \
        )

        ### Plot
        
        for ODF in ODF_plots:

            # Setup axis
            geo, prj = sfplt.getprojection(rotation=55+180, inclination=50)   
            W = 0.18 # axis width
            axin = plt.axes([ODF['axloc'][0],ODF['axloc'][1], W,W], projection=prj) #, transform=ax.transData)
            axin.set_global()
            axin.set_title(ODF['title'], ma='center', fontsize=FSANNO-0.5)
                    
            # Plot ODF
            sfplt.plotODF(ODF['xlm']*norm, lm, axin, cmap=ODF['cmap'], lvlset=(np.linspace(0.1,ODF['lvlmax'],7), lambda x,p:'%.1f'%x), showcb=False, nchunk=0)

    ### Limits, ticks, axis labels

    plt.sca(ax)
    
    plt.xlabel(r'$\hat{%s}_2^0$'%(axname))
    plt.ylabel(r'$\hat{%s}_4^0$'%(axname), labelpad=-2)

    plt.xlim(xlims)
    plt.ylim(ylims)

    ### Save figure

    fout = 'texplots/%s-%i.pdf'%(SELFNAME, ii)
    print('Saving %s'%(fout))
    plt.savefig(fout, transparent=True, dpi=250)
    plt.close()
    
os.system('cd texplots; pdflatex state-space-validation.tex; cd ..;')

