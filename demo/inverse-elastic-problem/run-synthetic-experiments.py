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
import math

from inverseproblem import * # Routines for solving inverse problems
from Cij import *            # Monocrystal elastic parameters

sys.path.insert(0, '..')
from specfabpy import specfabpy as sf 
lm, nlm_len = sf.init(4) # L=4 suffices here

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from plottools import * # For plotting ODFs etc.

#--------------------       
# Flags
#--------------------

VERBOSE_STATUS    = 1 # Verbose status 
VERBOSE_INVERSION = 0 # Verbose inversion 

#--------------------
# Synthetic experiments
#--------------------

x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

# Horizontal single maximum
a2_sm = np.tensordot(x,x, axes=0)
a4_sm = np.tensordot(a2_sm,a2_sm, axes=0) # <c^4> (aka. a^(4))
nlm_sm = sf.a4_to_nlm(a4_sm) 

# Horizontal double maximum
q = y 
q = 0.3*x + y; q /= np.linalg.norm(q) 
a2_smq = np.tensordot(q,q, axes=0)
a4_smq = np.tensordot(a2_smq,a2_smq, axes=0)
a4_dm = a4_sm/2 + a4_smq/2
nlm_dm = sf.a4_to_nlm(a4_dm)

# Vertical girdle
x2 = np.tensordot(x,x, axes=0); x4 = np.tensordot(x2,x2,axes=0)
z2 = np.tensordot(z,z, axes=0); z4 = np.tensordot(z2,z2,axes=0)
x2z2_2 = np.tensordot(x2+z2,x2+z2,axes=0)
xz, zx = np.tensordot(x,z, axes=0), np.tensordot(z,x, axes=0)
xz_zx_2 = np.tensordot(xz+zx,xz+zx,axes=0)
a4_gd = x4/4 + z4/4 + (x2z2_2)/8 + xz_zx_2/8
nlm_gd = sf.a4_to_nlm(a4_gd)

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
theta_xy, phi_xy = np.zeros(N),np.zeros(N)
theta_xz, phi_xz = np.zeros(N),np.zeros(N)

for nn, ang in enumerate(np.linspace(0,360,N)):
    Rz = R.from_euler('z', ang, degrees=True).as_matrix()
    Ry = R.from_euler('y', ang, degrees=True).as_matrix()
    v_xy = np.matmul(Rz,x)
    v_xz = np.matmul(Ry,x)
    theta_xy[nn], phi_xy[nn] = cart2sph(v_xy)
    theta_xz[nn], phi_xz[nn] = cart2sph(v_xz)


for ii, exprkey in enumerate(EXPERIMENTS):

    expr = EXPERIMENTS[exprkey]
    
    print('================')
    print(' EXPERIMENT: %s'%(expr['name']))
    print('================')

    #--------------------       
    # Modelled velocities from which to infer ODF
    #--------------------
    
    nlm = expr['nlm']
    vi_xy = get_vjmap(nlm, alpha,g_B68,rho, theta_xy,phi_xy) # forward modelled velocities
    vi_xz = get_vjmap(nlm, alpha,g_B68,rho, theta_xz,phi_xz) # forward modelled velocities

    print('*** Inferring ODF given V and g...')
    observations_xy = (vi_xy,theta_xy,phi_xy)
    observations_xz = (vi_xz,theta_xz,phi_xz)
    (nlm_xy_est, vi_xy_est, dvi_xy_est) = ip.infer_nlm(observations_xy, alpha, g_B68, beta, eta, use_angular_anomalies=False)
    (nlm_xz_est, vi_xz_est, dvi_xz_est) = ip.infer_nlm(observations_xz, alpha, g_B68, beta, eta, use_angular_anomalies=False)

    #--------------------
    # Plot results
    #--------------------

    ### S_2 projection

    inclination = 50 # view angle
    rot = -40 # view angle
    prj = ccrs.Orthographic(rot, 90-inclination)
    geo = ccrs.Geodetic()

    ### Setup figure, panels, etc.
    
    # Colors, lw, etc.
    c_xy  = '#1f78b4'
    cl_xy = c_xy # '#a6cee3'
    c_xz  = '#e31a1c'
    cl_xz = c_xz #'#fb9a99'

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
    
    axODF.set_global(); axODF_xy.set_global(); axODF_xz.set_global()

    ### ODFs
        
    tickintvl, lvls = 3, np.linspace(0.0,0.6,7)
    kwargs_ODF = {'cmap':'Greys', 'cbaspect':8.5, 'cbfrac':0.065,  'cborientation':'horizontal', 'lvls':lvls, 'tickintvl':tickintvl}
    plot_ODF(nlm, lm_L4, ax=axODF, cblabel='$\psi_{\mathrm{}}/N$', **kwargs_ODF)
#    qlatd, qlond = get_deg(qlat, qlon)
#    axODF_obs.plot(qlond, qlatd, ls='none', marker='x', markersize=0.6, c=c_caxes, transform=geo) 
    
    plot_ODF(nlm_xy_est, lm_L4, ax=axODF_xy, cblabel='$\psi_{\mathrm{}}/N$', **kwargs_ODF)
    plot_ODF(nlm_xz_est, lm_L4, ax=axODF_xz, cblabel='$\psi_{\mathrm{}}/N$', **kwargs_ODF)

    # Set ODF axes for reference
    plot_unitaxes(axODF,    geo)
    plot_unitaxes(axODF_xy, geo)
    plot_unitaxes(axODF_xz, geo)

    # Plot sampling directions
    kwargs = {'marker':'o', 'ms':2, 'transform':geo}
    for lat,lon in zip(theta_xy,phi_xy): 
        t,p = get_deg(lat, lon)
        axODF_xy.plot(p-180, t-90, color=c_xy, **kwargs)
    for lat,lon in zip(theta_xz,phi_xz): 
        t,p = get_deg(lat, lon)
        axODF_xz.plot(p-180, t-90, color=c_xz, **kwargs)

    # Titles
    pad = 10
    axODF.set_title(r'True ODF', pad=pad, fontsize=FS)
    axODF_xy.set_title('Inferred,\n $x$--$y$ plane', pad=pad, fontsize=FS)
    axODF_xz.set_title('Inferred,\n $x$--$z$ plane', pad=pad, fontsize=FS)

    # Panel no.
    kwargs_ODFpanelno = {'frameon':True, 'alpha':1.0, 'fontsize':FS}
    dx = +0.14 
    yloc=1.125 + 1*0.26
    writeSubplotLabel(axODF,    2, r'{\bf d}', bbox=(-0.46+dx,yloc), **kwargs_ODFpanelno)
    writeSubplotLabel(axODF_xy, 2, r'{\bf e}', bbox=(-0.31+dx,yloc), **kwargs_ODFpanelno)
    writeSubplotLabel(axODF_xz, 2, r'{\bf f}', bbox=(-0.28+dx,yloc), **kwargs_ODFpanelno)

    ### Velocity anomalies

    phideg = np.rad2deg(phi_xy)

    # Measured
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
    vi_xy_est = get_vjmap(nlm_xy_est, alpha, g_B68, rho, theta_xy, phi_xy) 
    vi_xz_est = get_vjmap(nlm_xz_est, alpha, g_B68, rho, theta_xz, phi_xz) 
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
    writeSubplotLabel(axP,  2, r'{\bf a}', frameon=True, alpha=1.0, fontsize=FS)
    writeSubplotLabel(axS1, 2, r'{\bf b}', frameon=True, alpha=1.0, fontsize=FS)
    writeSubplotLabel(axS2, 2, r'{\bf c}', frameon=True, alpha=1.0, fontsize=FS)

    # Axis labels
    axP.set_ylabel( r'$V_{\mathrm{qP}}$ ($\SI{}{\metre\per\second}$)')
    axS1.set_ylabel(r'$V_{\mathrm{qS1}}$ ($\SI{}{\metre\per\second}$)')
    axS2.set_ylabel(r'$V_{\mathrm{qS2}}$ ($\SI{}{\metre\per\second}$)')
    axS2.set_xlabel(r'$\phi$ or $\theta+\ang{90}$ ($\SI{}{\degree}$)')

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
     
