# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Script for running synthetic experiments: inferring the ODF from forward-modelled phase velocities *given* a prescribed ODF
These are the "additional experiments" in Rathmann et al. (2022).
"""

import copy, sys, code # code.interact(local=locals())

import numpy as np 
np.set_printoptions(edgeitems=30, linewidth=1000, formatter=dict(float=lambda x: "%+.6g" % x))
import numpy.linalg as linalg

from scipy.spatial.transform import Rotation as R

from inverseproblem import * # Routines for solving inverse problems
from Cij import *            # Monocrystal elastic parameters

from specfabpy import specfab as sf 
from specfabpy import discrete as sfdsc 
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

L = 4
lm, nlm_len = sf.init(L) # L=4 suffices here

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#--------------------       
# Flags
#--------------------

VERBOSE_STATUS    = 1 # Verbose status 
VERBOSE_INVERSION = 0 # Verbose inversion 

#--------------------
# Synthetic experiments
#--------------------

nlm_sm = sf.nlm_ideal([1,0,0], 0, L) # Horizontal single maximum
nlm_dm = (nlm_sm + sf.nlm_ideal([0,1,0], 0, L))/2 # Horizontal double maximum
nlm_gd = sf.rotate_nlm(sf.rotate_nlm(sf.nlm_ideal([0,0,1], np.pi/2, L), np.pi/2, 0), 0, np.pi/2) # Vertical girdle

EXPERIMENTS = {
    'singlemax': {'nlm':nlm_sm, 'name':'Horizontal single maximum'}, \
    'doublemax': {'nlm':nlm_dm, 'name':'Horizontal double maximum'}, \
    'girdle':    {'nlm':nlm_gd, 'name':'Vertical girdle'}, \
}

#--------------------       
# Monocrystal and homogenization parameters
#--------------------

# Monocrystal elastic parameters
g_B68 = Cij_to_g(Cij['Bennett1968']) # = (lam,mu,gam, Elam,Emu,Egam)
g_B68 = g_B68[[0,1,3,4,5]] # Routines below depend only on independent parameters (recall that gam=lam+2*mu)

# Reuss--Voigt weight of homogenization model
#alpha = 0 # Reuss
#alpha = 1 # Voigt
alpha = 0.5 # Hill average

rho = 917 # Mass density

#--------------------       
# Init
#--------------------

ip = InverseProblem(rho=rho, verbose=VERBOSE_INVERSION)

# Misfit weights (qP, qS1, qS2)
beta = [1,1,1] # => equal uncertainty on P- ans S-wave velocities

# Lagrange multipliers for imposing inequality constraints a_1>0 and \tilde{a}_1>0
#eta = [0, 2.5e3] # default eigenvalue regularization 
eta = [0, 0] # no eigenvalue regularization

#--------------------
# Infer ODFs 
#--------------------

N = 20 # number of equi-spaced sampling directions
v_xy, colat_xy, lon_xy = np.zeros((N,3)), np.zeros(N), np.zeros(N)
v_xz, colat_xz, lon_xz = np.zeros((N,3)), np.zeros(N), np.zeros(N)

x,y,z = np.eye(3)

for nn, ang in enumerate(np.linspace(0,360,N)):
    Rz = R.from_euler('z', ang, degrees=True).as_matrix()
    Ry = R.from_euler('y', ang, degrees=True).as_matrix()
    v_xy[nn,:] = np.matmul(Rz,x)
    v_xz[nn,:] = np.matmul(Ry,x)
    _, colat_xy[nn], lon_xy[nn] = cart2sph_wrap(v_xy[nn,:]) # colat, lon
    _, colat_xz[nn], lon_xz[nn] = cart2sph_wrap(v_xz[nn,:]) # colat, lon


for ii, exprkey in enumerate(EXPERIMENTS):

    expr = EXPERIMENTS[exprkey]
    
    print('================')
    print(' EXPERIMENT: %s'%(expr['name']))
    print('================')

    #--------------------       
    # Modelled velocities from which to infer ODF
    #--------------------
    
    nlm = expr['nlm']
    vi_xy = get_vi_map(nlm, alpha,g_B68,rho, colat_xy,lon_xy) # forward modelled velocities
    vi_xz = get_vi_map(nlm, alpha,g_B68,rho, colat_xz,lon_xz) # forward modelled velocities

    print('*** Inferring ODF given V and g...')
    observations_xy = (vi_xy,colat_xy,lon_xy)
    observations_xz = (vi_xz,colat_xz,lon_xz)
    (nlm_xy_est, vi_xy_est, dvi_xy_est) = ip.infer_nlm(observations_xy, alpha, g_B68, beta, eta, use_angular_anomalies=False)
    (nlm_xz_est, vi_xz_est, dvi_xz_est) = ip.infer_nlm(observations_xz, alpha, g_B68, beta, eta, use_angular_anomalies=False)

    #--------------------
    # Plot results
    #--------------------

    ### S_2 projection

    geo, prj = sfplt.getprojection(rotation=50, inclination=50)

    ### Setup figure, panels, etc.
    
    # Colors, lw, etc.
    cl_xy = c_xy = sfplt.c_blue
    cl_xz = c_xz = sfplt.c_red

    scale = 0.5
    fig = plt.figure(figsize=(14*scale,14.5*scale))

    gs_master = gridspec.GridSpec(2, 1, height_ratios=[1,0.3]) 
    a=0.11
    gs_master.update(top=0.99, bottom=0.07, left=a, right=1-a*0.25, hspace=0.4)

    gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_master[0,0], hspace=0.25)
    axP  = fig.add_subplot(gs[0,:])
    axS1 = fig.add_subplot(gs[1,:], sharex=axP)
    axS2 = fig.add_subplot(gs[2,:], sharex=axP)

    gs = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs_master[1,0],  width_ratios=[-0.55,1,1,1,1.65])
    ii0=1
    axODF    = fig.add_subplot(gs[0, ii0+0], projection=prj)
    axODF_xy = fig.add_subplot(gs[0, ii0+1], projection=prj)
    axODF_xz = fig.add_subplot(gs[0, ii0+2], projection=prj)
    
    axODF.set_global()
    axODF_xy.set_global()
    axODF_xz.set_global()

    ### ODFs
        
    kwargs = dict(lvlset=(np.linspace(0.0,0.6,7), lambda x,p:'%.1f'%x), cbtickintvl=3, cbaspect=8.5, cbfraction=0.065)
    sfplt.plotODF(nlm, lm_L4, axODF, cblabel='$\psi_{\mathrm{}}/N$', **kwargs)
    sfplt.plotODF(nlm_xy_est, lm_L4, axODF_xy, cblabel='$\psi_{\mathrm{}}/N$', **kwargs)
    sfplt.plotODF(nlm_xz_est, lm_L4, axODF_xz, cblabel='$\psi_{\mathrm{}}/N$', **kwargs)

    # Set ODF axes for reference
    for axi in (axODF, axODF_xy, axODF_xz): 
        sfplt.plotcoordaxes(axi, geo, axislabels='vuxi', color=sfplt.c_dred)  

    # Plot sampling directions
    kwargs = {'marker':'o', 'ms':2, 'transform':geo}
    for nn in range(N): 
        sfplt.plotS2point(axODF_xy, v_xy[nn,:], color=c_xy, **kwargs)
        sfplt.plotS2point(axODF_xz, v_xz[nn,:], color=c_xz, **kwargs)

    # Titles
    pad = 10
    axODF.set_title(r'True ODF', pad=pad, fontsize=FS)
    axODF_xy.set_title('Inferred,\n $x$--$y$ plane', pad=pad, fontsize=FS)
    axODF_xz.set_title('Inferred,\n $x$--$z$ plane', pad=pad, fontsize=FS)

    # Panel no.
    kwargs_ODFpanelno = {'frameon':True, 'alpha':1.0, 'fontsize':FS}
    dx = +0.14 
    yloc=1.125 + 1*0.26
    sfplt.panellabel(axODF,    2, r'{\bf d}', bbox=(-0.46+dx,yloc), **kwargs_ODFpanelno)
    sfplt.panellabel(axODF_xy, 2, r'{\bf e}', bbox=(-0.31+dx,yloc), **kwargs_ODFpanelno)
    sfplt.panellabel(axODF_xz, 2, r'{\bf f}', bbox=(-0.28+dx,yloc), **kwargs_ODFpanelno)

    ### Velocity anomalies

    phideg = np.rad2deg(lon_xy)

    # True values
    (vP_xy, vS1_xy, vS2_xy) = vi_xy
    (vP_xz, vS1_xz, vS2_xz) = vi_xz
    kwargs = {'marker':'o', 'ls':'none', 'ms':5.5}
    axP.plot( phideg, vP_xy,  c=cl_xy, label=r'True ODF, $x$--$y$ plane', **kwargs)
    axP.plot( phideg, vP_xz,  c=cl_xz, label=r'True ODF, $x$--$z$ plane', **kwargs)
    axS1.plot(phideg, vS1_xy, c=cl_xy, **kwargs)
    axS1.plot(phideg, vS1_xz, c=cl_xz, **kwargs)
    axS2.plot(phideg, vS2_xy, c=cl_xy, **kwargs)
    axS2.plot(phideg, vS2_xz, c=cl_xz, **kwargs)

    # Predicted velocities using observed ODF and g_B68
    vi_xy_est = get_vi_map(nlm_xy_est, alpha, g_B68, rho, colat_xy, lon_xy) 
    vi_xz_est = get_vi_map(nlm_xz_est, alpha, g_B68, rho, colat_xz, lon_xz) 
    axP.plot( phideg, vi_xy_est[0], '--', c=c_xy, label=r'Inferred ODF, $x$--$y$ plane')
    axP.plot( phideg, vi_xz_est[0], '--', c=c_xz, label=r'Inferred ODF, $x$--$z$ plane')
    axS1.plot(phideg, vi_xy_est[1], '--', c=c_xy)
    axS1.plot(phideg, vi_xz_est[1], '--', c=c_xz)
    axS2.plot(phideg, vi_xy_est[2], '--', c=c_xy) 
    axS2.plot(phideg, vi_xz_est[2], '--', c=c_xz)

    # Legend
    legkwargs = {'frameon':False, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.95, 'ncol':1, 'handlelength':1.2, 'columnspacing':1, 'labelspacing':0.3}
    hleg = axP.legend(loc=1, fontsize=FS-0.5, bbox_to_anchor=(0.97,-2.94), **legkwargs)
    hleg.get_frame().set_linewidth(0.8);

    # Panel no.
    sfplt.panellabel(axP,  2, r'{\bf a}', frameon=True, alpha=1.0, fontsize=FS)
    sfplt.panellabel(axS1, 2, r'{\bf b}', frameon=True, alpha=1.0, fontsize=FS)
    sfplt.panellabel(axS2, 2, r'{\bf c}', frameon=True, alpha=1.0, fontsize=FS)

    # Axis labels
    axP.set_ylabel( r'$V_{\mathrm{qP}}$ ($\SI{}{\metre\per\second}$)')
    axS1.set_ylabel(r'$V_{\mathrm{qS1}}$ ($\SI{}{\metre\per\second}$)')
    axS2.set_ylabel(r'$V_{\mathrm{qS2}}$ ($\SI{}{\metre\per\second}$)')
    axS2.set_xlabel(r'$\phi$, $\theta$ ($\SI{}{\degree}$)')

    # y-axis limits
    yticks = np.arange(-200,200+1e-5,50)

    # Common x-axis
    xticks = np.arange(0,360+1e-5,45)
    axS2.set_xticks(xticks[::2])
    axS2.set_xticks(xticks[::], minor=True)
    axS2.set_xlim([0,360])
    plt.setp(axP.get_xticklabels(), visible=False) 
    plt.setp(axS1.get_xticklabels(), visible=False) 

    ### Save figure
    
    fname = 'plots/synthetic-%s.png'%(exprkey)
    print('*** Saving inversion results to %s'%(fname))
    plt.savefig(fname, dpi=150)
    plt.close()
     
