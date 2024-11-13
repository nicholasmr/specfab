# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Main script for producing the plots and tables published in Rathmann et al. (2022), doi:10.1098/rspa.2022.0574  

Results may vary slightly depending on numpy and scipy versions.
"""

import copy, sys, code # code.interact(local=locals())
import numpy as np 
np.set_printoptions(edgeitems=30, linewidth=1000, formatter=dict(float=lambda x: "%+.6g" % x))
import numpy.linalg as linalg

from experiments import *    # Routines for handling Lutz et al. (2022) experiments and more
from inverseproblem import * # Routines for solving inverse problems
from Cij import *            # Monocrystal/grain elastic parameters

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#--------------------       
# Flags
#--------------------

VERBOSE_STATUS    = 1 # Verbose status 
VERBOSE_INVERSION = 0 # Verbose inversion 

MAKE_GAW_STATS = 0

#--------------------
# Lutz et al. (2022) experiments to process
#--------------------

EXPERIMENTS = [3,7,10]
#EXPERIMENTS = [7,] # debug

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
# Inversion parameters
#--------------------

# Misfit weights (qP, qS1, qS2)
beta = [1,1,1]
#beta = [1,0,0] # Consider only P-wave data (just for testing its effect...)

# Lagrange multipliers for imposing inequality constraints a_1>0 and \tilde{a}_1>0
eta = {\
    '003': [0, 1.4e3], \
    '007': [0, 2.5e3], \
    '010': [0, 0.5e3], \
}

g0 = [g_B68[0],g_B68[1], 1,1,1] # initial guess for estimating elastic parameter vector "g"

#--------------------       
# Init
#--------------------

lutz = Lutz_etal_2022(rho=rho, verbose=VERBOSE_INVERSION)
ip   = InverseProblem(rho=rho, verbose=VERBOSE_INVERSION)

# For later printing latex tables
eigvals   = np.zeros((3+2,3,3)) # Eigenvalues for Lutz et al. samples
g_RVH     = np.ones((5,3,3))   # Inferred elastic parameters for Lutz et al. samples using Reuss(R)/Voigt(V)/Hill(H) homogenizations
g_RVH_gaw = np.ones((5,3,3))   # Same, but for grain-area-weighted (gaw) ODFs

#--------------------
# Process experiments
#--------------------

for ii, exprnum in enumerate(EXPERIMENTS):

    print('================')
    print(' EXPERIMENT %03i'%(exprnum))
    print('================')

    #--------------------       
    # Load observed ODF and velocities
    #--------------------

    (nlm_obs, qlat,qlon) = lutz.get_nlm(exprnum) 
    (theta,phi, vi_obs, visig_obs) = lutz.get_velocities(exprnum)
    (vP_obs,vS1_obs,vS2_obs) = vi_obs # unpack
    (vP_obs_sig,vS1_obs_sig,vS2_obs_sig) = visig_obs # unpack
    
    #--------------------
    # Best-fit elastic parameters, g_R22, given the true ODF and phase-velocity observations
    #--------------------

    observations = (vi_obs,nlm_obs, theta,phi)

    print('*** Inferring g_R22 given obs. ODF, alpha=%.2f'%(alpha))    
    (g_R22, vi_p, dvi_p) = ip.infer_g(observations, alpha,beta, g0)
    print('----------------------')
    print('g_B68 = ', g_B68)
    print('g_R22 = ', g_R22)
#    print('rel. pct. = ', np.divide(g_R22-g_B68,g_B68)*100)
    print('----------------------')
    
    print('*** Inferring g_R22 given obs. ODF for alpha=0,1,1/2 (for latex table)')    
    # Voigt, Reuss, and Hill parameters for constructing latex table
    (g_RVH[:,ii,0], _,_) = ip.infer_g(observations, 0.0, beta, g0)
    (g_RVH[:,ii,1], _,_) = ip.infer_g(observations, 1.0, beta, g0)
    (g_RVH[:,ii,2], _,_) = ip.infer_g(observations, 0.5, beta, g0)

    ### Grain area weighted versions (gaw)
    if MAKE_GAW_STATS:
        print('*** Inferring g_R22 given obs. ODF with area-weigted c-axes for alpha=0,1,1/2 (for latex table)')    
        (nlm_obs_gaw, _,_) = lutz.get_nlm(exprnum, weighted=True)
        observations = (vi_obs,nlm_obs_gaw, theta,phi)
        (g_RVH_gaw[:,ii,0], _,_) = ip.infer_g(observations, 0.0, beta, g0)
        (g_RVH_gaw[:,ii,1], _,_) = ip.infer_g(observations, 1.0, beta, g0)
        (g_RVH_gaw[:,ii,2], _,_) = ip.infer_g(observations, 0.5, beta, g0) 

    #--------------------       
    # Infer ODF given g_B68
    #--------------------
    
    observations = (vi_obs, theta,phi)
    eta_ = eta['%03i'%(exprnum)]

    print('*** Inferring ODF given g_B68...')
    (nlm_B68, vi_B68, dvi_B68) = ip.infer_nlm(observations, alpha, g_B68, beta, eta_, use_angular_anomalies=True)

    print('*** Inferring ODF given g_R22...')
    (nlm_R22, vi_R22, dvi_R22) = ip.infer_nlm(observations, alpha, g_R22, beta, eta_, use_angular_anomalies=False)

    #--------------------
    # Print eigen values
    #--------------------

    if VERBOSE_STATUS: print('*** Eigenvalues are')
    print('----------------------')
    names = [r'Observed', r'g_B68', r'g_R22']
    for nn, nlm_ in enumerate((nlm_obs, nlm_B68, nlm_R22)):
        lami = np.sort(linalg.eigvals(sf.a2(nlm_)))
        _,_,_,_,_,_, Lami = sf.a4_eigentensors(nlm_)
        eigvals[:,ii,nn] = np.concatenate((lami,[np.amin(Lami),np.amax(Lami)]))
        if VERBOSE_STATUS: print(r'%8s // a2 eigvals: %.2f, %.2f, %.2f (sum %.2f) // a4 eigvals: %.2f, %.2f, %.2f, %.2f, %.2f, %.2f (sum %.2f)'%(names[nn], lami[0], lami[1], lami[2], lami[0]+lami[1]+lami[2], \
                                                                                                                                                            Lami[0],Lami[1],Lami[2],Lami[3],Lami[4],Lami[5], np.sum(Lami)))
    print('----------------------')
    
    #--------------------
    # Plot results
    #--------------------

    ### S_2 projection

    geo, prj = sfplt.getprojection(rotation=50-90, inclination=50)

    ### Setup figure, panels, etc.
    
    # Colors, lw, etc.
    c_obs   = 'k'
    c_B68   = '#fb9a99'
    c_R22   = '#e31a1c'
    c_p     = '#1f78b4'
    c_pB68  = '#a6cee3'
    c_caxes = '#33a02c'

    scale = 0.5
    fig = plt.figure(figsize=(14*scale,14.5*scale))

    gs_master = gridspec.GridSpec(2, 1, height_ratios=[1,0.3]) 
    a=0.11
    gs_master.update(top=0.975, bottom=0.07, left=a, right=1-a*0.25, hspace=0.4)

    gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_master[0,0], hspace=0.25)
    axP  = fig.add_subplot(gs[0,:])
    axS1 = fig.add_subplot(gs[1,:], sharex=axP)
    axS2 = fig.add_subplot(gs[2,:], sharex=axP)

    gs = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs_master[1,0],  width_ratios=[-0.55,1,1,1,1.65])
    ii0=1
    axODF_obs = fig.add_subplot(gs[0, ii0+0], projection=prj)
    axODF_B68 = fig.add_subplot(gs[0, ii0+1], projection=prj)
    axODF_R22 = fig.add_subplot(gs[0, ii0+2], projection=prj)
    
    axODF_obs.set_global()
    axODF_B68.set_global()
    axODF_R22.set_global()

    ### ODFs
        
    kwargs = dict(lvlset=(np.linspace(0.0,0.6,7), lambda x,p:'%.1f'%x), cbtickintvl=3, cbaspect=8.5, cbfraction=0.065)
    sfplt.plotODF(nlm_obs, lm_L4, axODF_obs, cblabel='$\psi_{\mathrm{}}/N$', **kwargs)
    sfplt.plotODF(nlm_B68, lm_L4, axODF_B68, cblabel='$\psi_{\mathrm{}}/N$', **kwargs)
    sfplt.plotODF(nlm_R22, lm_L4, axODF_R22, cblabel='$\psi_{\mathrm{}}/N$', **kwargs)   

    axODF_obs.plot(np.rad2deg(qlon), np.rad2deg(qlat), ls='none', marker='x', markersize=0.6, c=c_caxes, transform=geo) 

    # Set ODF axes for reference
    for axi in (axODF_obs, axODF_B68, axODF_R22): 
        sfplt.plotcoordaxes(axi, geo, axislabels='vuxi', color=sfplt.c_dred)  
        
    # Titles
    pad = 10
    axODF_obs.set_title(r'observed, $\hat{\vb*{\psi}}_{\mathrm{obs}}$', pad=pad, fontsize=FS)
    axODF_B68.set_title('$\\hat{\\vb*{\psi}}$ inferred\n given $\\vb{g}_{\\mathrm{B68}}$', pad=pad, fontsize=FS)
    axODF_R22.set_title('$\\hat{\\vb*{\psi}}$ inferred\n given $\\vb{g}_{\\mathrm{opt}}$', pad=pad, fontsize=FS)

    # Panel no.
    kwargs_ODFpanelno = {'frameon':False, 'alpha':1.0, 'fontsize':FS}
    dx = +0.14 
    yloc=1.125 + 2*0.26
    sfplt.panellabel(axODF_obs, 2, r'(\textit{d})', bbox=(-0.35+dx,yloc), **kwargs_ODFpanelno)
    sfplt.panellabel(axODF_B68, 2, r'(\textit{e})', bbox=(-0.35+dx,yloc), **kwargs_ODFpanelno)
    sfplt.panellabel(axODF_R22, 2, r'(\textit{f})', bbox=(-0.35+dx,yloc), **kwargs_ODFpanelno)

    ### Velocity anomalies

    phideg = np.rad2deg(phi)

    # Measured
    errbarkwargs = {'fmt':'o', 'c':c_obs, 'capsize':1.5, 'ms':4.5, 'lw':1, 'zorder':0}
    axP.errorbar( phideg, vP_obs,  yerr=vP_obs_sig,  label=r'observed', **errbarkwargs)
    axS1.errorbar(phideg, vS1_obs, yerr=vS1_obs_sig, **errbarkwargs)
    axS2.errorbar(phideg, vS2_obs, yerr=vS2_obs_sig, **errbarkwargs)

    # Predicted velocities using observed ODF and g_B68
    vi_pB68 = get_vi_map(nlm_obs, alpha, g_B68, rho, theta, phi) 
    axP.plot( phideg, vi_pB68[0], '--', c=c_pB68, label=r'prediction given $\vb{g}_{\mathrm{B68}}$ and $\hat{\vb*{\psi}}_{\mathrm{obs}}$')
    axS1.plot(phideg, vi_pB68[1], '--', c=c_pB68)
    axS2.plot(phideg, vi_pB68[2], '--', c=c_pB68) 

    # g_opt | ODF_obs
    axP.plot( phideg, vi_p[0], '-', c=c_p, label=r'$\vb{g}_{\mathrm{opt}}$ inferred given $\hat{\vb*{\psi}}_{\mathrm{obs}}$')
    axS1.plot(phideg, vi_p[1], '-', c=c_p)
    axS2.plot(phideg, vi_p[2], '-', c=c_p)

    # ODF | g_B68
    axP.plot( phideg, vi_B68[0], '--', c=c_B68, label=r'$\hat{\vb*{\psi}}$ inferred given $\vb{g}_{\mathrm{B68}}$')
    axS1.plot(phideg, vi_B68[1], '--', c=c_B68)
    axS2.plot(phideg, vi_B68[2], '--', c=c_B68)

    # ODF | g_R22
    axP.plot( phideg, vi_R22[0], '-', c=c_R22, label=r'$\hat{\vb*{\psi}}$ inferred given $\vb{g}_{\mathrm{opt}}$')
    axS1.plot(phideg, vi_R22[1], '-', c=c_R22)
    axS2.plot(phideg, vi_R22[2], '-', c=c_R22)

    # Legend
    legkwargs = {'frameon':False, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.95, 'ncol':1, 'handlelength':1.2, 'columnspacing':1, 'labelspacing':0.3}
    handles, labels = axP.get_legend_handles_labels()
#    order = [0,1,4,2,3]
    order = [4,0,1,2,3]
    hleg = axP.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc=1, fontsize=FS-0.5, bbox_to_anchor=(1.05,-2.94), **legkwargs)
    hleg.get_frame().set_linewidth(0.8);

    # Panel no.
    bbox=(-0.14,1.23)
    sfplt.panellabel(axP,  2, r'(\textit{a})', frameon=False, alpha=1.0, fontsize=FS, bbox=bbox)
    sfplt.panellabel(axS1, 2, r'(\textit{b})', frameon=False, alpha=1.0, fontsize=FS, bbox=bbox)
    sfplt.panellabel(axS2, 2, r'(\textit{c})', frameon=False, alpha=1.0, fontsize=FS, bbox=bbox)

    # Axis labels
    axP.set_ylabel( r'$V_{\mathrm{qP}}$ ($\SI{}{\metre\per\second}$)')
    axS1.set_ylabel(r'$V_{\mathrm{qS1}}$ ($\SI{}{\metre\per\second}$)')
    axS2.set_ylabel(r'$V_{\mathrm{qS2}}$ ($\SI{}{\metre\per\second}$)')
    axS2.set_xlabel(r'$\phi$ ($\SI{}{\degree}$)')

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
    
    fsuffix = 'Lutz%03i'%(exprnum)
    fname = 'plots/infer-%s-alpha=%.2f.png'%(fsuffix, alpha)
#    fname = 'plots/infer-%s-alpha=%.2f.pdf'%(fsuffix, alpha)
    print('*** Saving inversion results to %s'%(fname))
    plt.savefig(fname, dpi=300)
    
#--------------------
# Python exports
#--------------------
    
### Inferred "g" vectors
print('*** Best-fit g vector components are:')
print('-------------------------------------')
for nn,exprnum in enumerate([3,7,10]): 
    for ii,alphastr in enumerate(['0','1','5']):
        g_ = g_RVH[:,nn,ii]
        print("\t'%03i_%s': [%.4e, %.4e, %.4f, %.4f, %.4f],\\"%(exprnum,alphastr, g_[0]*1e-9,g_[1]*1e-9,g_[2],g_[3],g_[4]))
print('-------------------------------------')

#--------------------
# Latex exports
#--------------------
    
### Eigenvalues
fout = "latex/eigenvalues.tex"
print('*** Saving eigenvalues to %s'%(fout))
ftable = open(fout, "w")
names = [r'Observed', r'Inferred, $\vb{g}_{\mathrm{B68}}$', r'Inferred, $\vb{g}_{\mathrm{opt}}$']
for nn in np.arange(3): # ODF 
    ftable.write('%33s'%(names[nn]))
    for ii in np.arange(3): # sample
        ftable.write(r'& %.2f & %.2f & %.2f & %.2f  '%(eigvals[0,ii,nn], eigvals[1,ii,nn], eigvals[2,ii,nn], eigvals[3,ii,nn]))
    ftable.write('\\\\ \n')
ftable.close()


### Monocrystal parameters for **grain-area-weigted** ODF
samplenames = [r'Sample 003', r'Sample 007', r'Sample 010']
if MAKE_GAW_STATS:
    fout = "latex/crystalparams_gaw.tex"
    print('*** Saving variations in inferred monocrystal paramters for GAW ODFs to %s'%(fout))
    ftable = open(fout, "w")
    for nn in np.arange(3): # ODF 
        ftable.write('\\addlinespace[6pt]\n')
        ftable.write('\\multicolumn{1}{l}{\\bf{%s}} & \\\\ \n'%(samplenames[nn]))
        for ii,alphastr in enumerate(['0','1','1/2']):
            gamma     = g_RVH[0,nn,ii] + 2*g_RVH[1,nn,ii]
            gamma_gaw = g_RVH_gaw[0,nn,ii] + 2*g_RVH_gaw[1,nn,ii]
            dlam  = 100*(g_RVH_gaw[0,nn,ii]-g_RVH[0,nn,ii])/g_RVH[0,nn,ii]
            dmu   = 100*(g_RVH_gaw[0,nn,ii]-g_RVH[0,nn,ii])/g_RVH[0,nn,ii]
            dgam  = 100*(gamma_gaw-gamma)/gamma
            dElam = 100*(g_RVH_gaw[0,nn,ii]-g_RVH[0,nn,ii])/g_RVH[0,nn,ii]
            dEmu  = 100*(g_RVH_gaw[0,nn,ii]-g_RVH[0,nn,ii])/g_RVH[0,nn,ii]
            dEgam = 100*(g_RVH_gaw[0,nn,ii]-g_RVH[0,nn,ii])/g_RVH[0,nn,ii]
            ftable.write('$\\alpha=%3s$ & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% \\\\ \n'%(alphastr, dlam, dmu, dgam, dElam, dEmu, dEgam))
    ftable.close()

# Largest difference in inferred parameters among the homogenization schemes
fout = "latex/crystalparams_maxvariance.tex"
print('*** Saving maximum variations found in g across inversions to %s'%(fout))
ftable = open(fout, "w")
# g_RVH = (param, lutz sample, R/V/H approx)
g_delta = [np.amax([  (np.amax(g_RVH[jj,ii,:])-np.amin(g_RVH[jj,ii,:]))/np.amin(g_RVH[jj,ii,:])  for ii in np.arange(3)]) for jj in np.arange(5)]
ftable.write('$\delta\lambda = \SI{%.1f}{\percent}$,\n' %(100*np.abs( g_delta[0] ))) 
ftable.write('$\delta\mu     = \SI{%.1f}{\percent}$,\n' %(100*np.abs( g_delta[1] ))) 
ftable.write('$\delta\Elam   = \SI{%.1f}{\percent}$,\n' %(100*np.abs( g_delta[2] ))) 
ftable.write('$\delta\Emu    = \SI{%.1f}{\percent}$,\n' %(100*np.abs( g_delta[3] ))) 
ftable.write('$\delta\Egam   = \SI{%.1f}{\percent}$.\n' %(100*np.abs( g_delta[4] )))
ftable.close()
 
