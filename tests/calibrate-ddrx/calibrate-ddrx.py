# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from scipy.interpolate import interp1d
import pickle
from progress.bar import Bar

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf

yr2s = 31556926
s2yr = 3.16887646e-8
c2k  = lambda degc: degc+273.15 # deg. C to deg. K

#---------------------
# Numerics
#---------------------

L_list = [4,6,8,20] # Spectral truncation used for calibration
L_list = [8,]

Nt = 200 # Number of integration steps

WITH_CDRX = 1

for L in L_list:

    #---------------------
    # Gamma_0 
    #---------------------

    Q = 1e5 # activation energy

    if L == 4:
        A = 5.0e20 
    if L == 6:
        A = 5.0e20 
    if L == 8:
        A = 4.5e20 
    if L == 20:
        A = 2.0*5.0e20 

    print('Gamma(T=-10)/Gamma(T=-20) = %e'%(sf.Gamma0(np.eye(3),c2k(-10),A,Q)/sf.Gamma0(np.eye(3),c2k(-20),A,Q)))

    #---------------------
    # Load EDC profile
    #---------------------

    PROFILE = 'EDC'

    if PROFILE == 'EDC':
        file = open("DOME_C.dump",'rb')
        (fab,temp,aux) = pickle.load(file)
        n20_0 = 0.14 # intial fabric
        strain_zz_stop=-0.96 # stop early, fabric only measure at EDC until roughly this point.

    # Geometry
    z_bed = -aux['H'] # depth of bed
    z0 = abs(fab['z'][0]) 
    z1 = aux['H']
    H = z1-z0 # ice column height as seen by parcel model
    print('z_bed, H = %f, %f'%(z_bed, H))

    # Fabric
    lami_meas, z_meas = fab['eig'], fab['z']

    # Temperature interpolant
    T, zT = temp['T'], temp['z']
    f_T = interp1d(zT,T)

    # Accumulation rate
    b = aux['b']

    #---------------------
    # Pure shear deformation 
    #---------------------

    class PureShear():
        def __init__(self, t_e, r=0, ax='z'): 
            self.t_e = float(t_e) # e-folding time scale for for parcel height reduction
            if ax=='z': self.Fpow = [(1+r)/2., (1-r)/2., -1] # r=0 => unconfined
        def lam(self, t):    return np.exp(t/self.t_e) # lambda(t)
        def F(self, t):      return np.diag(np.power(self.lam(t),self.Fpow)) # Deformation tensor
        def strain(self, t): return 0.5*( self.F(t) + np.transpose(self.F(t)) ) - np.diag([1,1,1]) # Strain tensor
        def strainzz2time(self,strain_zz): return -self.t_e*np.log(strain_zz+1) # time it takes to reach "strain_zz" strain with timescale t_e.
        def D(self): return 1/self.t_e * np.diag(self.Fpow) # Strain rate. Note that F is constructed such that W and eps are time-independant.
        def W(self): return np.diag([0,0,0]) # Spin

    #---------------------
    # Ice parcel model
    #---------------------
        
    # Assumes depth-constant (unconfined) vertical compression: the Nye's classical dome model.
    b /= yr2s
    t_e = 1/(b/H) # e-folding time for uniaxial compression (strainrate_zz = MeanAccum/H)
    ps = PureShear(t_e, r=0.0, ax='z') # r=0 => unconfined
        
    # Constants
    tend    = ps.strainzz2time(strain_zz_stop)
    timevec = np.linspace(0,tend,Nt)
    dt      = timevec[1]-timevec[0] # given t_e, if we want Nt time steps, this is the step size needed
    strain_zz = np.array([ps.strain(t)[-1,-1] for ii, t in enumerate(timevec)]) # for the constant time-step size, these are the vertical parcel strains as a function of time
    z_sf = H*(strain_zz+1) - H # Vertical parcel strains correspond to these depths
    z_sf += -z0 # effective offset of parcel model

    # Initialize arrays
    lm, nlm_len = sf.init(L)
    nlm = np.zeros((Nt, nlm_len), dtype=np.complex64) # array of expansion coefficients
    nlm[0,0] = 1/np.sqrt(4*np.pi) # normalized distribution
    nlm[0,3] = n20_0 # initial single max strength (fitted to match eigen value profile)
    lami_sf = np.zeros((Nt, 3)) # array of modelled eigen values
    lami_sf[0,:] = np.diag(sf.a2(nlm[0,:])) # a2 already diagional (in eigenbasis) for this mode of deformation, so simply extract the diagonal values
    gam, lam = np.zeros((Nt)), np.zeros((Nt))
        
    # Euler integration of fabric model
    D, W = ps.D(), ps.W() # strain-rate and spin tensors for mode of deformation
    T = f_T(z_sf)
    gam[0] = sf.Gamma0(D, c2k(T[0]), A, Q)
    def f_lam(D,T):
        m_lam, b_lam = 1.26e-3, 0.206
        eps_E = np.sqrt(0.5*np.einsum('ij,ji',D,D)) # sqrt(0.5 * D:D)
        return eps_E * (m_lam*T + b_lam)
    lam[0] = f_lam(D,T[0])
    
    with Bar('dt=%.3fyr, Nt=%i :: L=%i (nlm_len=%i) ::'%(dt*s2yr,Nt,L,nlm_len), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds') as bar:
        for nn in np.arange(1,Nt):
        
            nlm_prev = nlm[nn-1,:]
            
            M_LATROT = sf.dndt_LATROT(nlm_prev, D, W) # Lattice rotation operator (nlm_len x nlm_len matrix)

            S = D.copy() # S remains coaxial with D for this mode of deformation, so we need not calculate S from the bulk flow law (which in turn requires calculating the enhancement factors from the modelled fabric)
            gam[nn] = sf.Gamma0(D, c2k(T[nn]), A, Q) # DDRX decary rate magnitude
            M_DDRX = gam[nn] * sf.dndt_DDRX(nlm_prev, S) # DDRX operator (nlm_len x nlm_len matrix)

            if WITH_CDRX:
                lam[nn] = f_lam(D,T[nn]) 
                M_CDRX = lam[nn]*sf.dndt_CDRX(nlm_prev) # CDRX operator (nlm_len x nlm_len matrix)
            else:
                M_CDRX = 0*M_DDRX

            M_REG = sf.dndt_REG(nlm_prev, D) # Regularization operator (nlm_len x nlm_len matrix)

            M = M_LATROT + M_DDRX + M_CDRX + M_REG # Total fabric evolution operator (matrix)
            nlm[nn,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # Forward Euler step
            
            nlm[nn,:] = sf.apply_bounds(nlm[nn,:])    # Apply spectral bounds if needed 
            
            lami_sf[nn,:] = np.diag(sf.a2(nlm[nn,:])) # a2 already diagional (in eigen basis) for this mode of deformation, so simply extract the diagonal values
            
            bar.next()
            
    ### Estimate numerical accuracy 
    diff_abs = np.abs(100*(nlm[:,0]-nlm[0,0])/nlm[0,0])
    print(r'Numerical error: max(n00(t) - n00(t=0)) [%] = ', np.amax(diff_abs))

    #--------------
    # Plot results
    #--------------

    import matplotlib.pyplot as plt

    cred    = '#e31a1c'
    credl   = '#fb9a99'
    cgreen  = '#33a02c'
    cgreenl = '#b2df8a'
    cblue   = '#1f78b4'
    cbluel  = '#a6cee3'

    scale = 0.55
    fig = plt.figure(figsize=(14*scale,9*scale))
    gs = fig.add_gridspec(ncols=4, nrows=3, width_ratios=[1.8,1.3,1,1], right=0.97, top=0.95, hspace=0.55)
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[:, 2], sharey=ax1)
    ax3 = fig.add_subplot(gs[:, 3], sharey=ax1)

    ### Eigen vals
    ax1.scatter(lami_meas[:,0], z_meas, c=cgreenl)
    ax1.scatter(lami_meas[:,1], z_meas, c=cbluel)
    ax1.scatter(lami_meas[:,2], z_meas, c=credl)
    if 0:
        lw, ls = 1.5, '--'
        ax1.plot(lami[:,0],z[:], c=cgreen, ls=ls, lw=lw)
        ax1.plot(lami[:,1],z[:], c=cblue, ls=ls, lw=lw)
        ax1.plot(lami[:,2],z[:], c=cred, ls=ls, lw=lw)

    ax1.plot(lami_sf[:,0],z_sf, c=cgreen, ls='-', lw=2, clip_on=False)
    ax1.plot(lami_sf[:,1],z_sf, c=cblue, ls='--', lw=2, clip_on=False)
    ax1.plot(lami_sf[:,2],z_sf, c=cred, ls='-', lw=2, clip_on=False)

    ax1.set_xlabel(r'$\lambda_i$')
    ax1.set_xlim([0,1])

    ax1.set_ylabel(r'$z$ (m)')
    ax1.set_ylim([z_bed,0])
    ax1.set_yticks(-np.flipud(np.arange(0,H,500)))
    ax1.set_yticks(-np.flipud(np.arange(0,H,250)), minor=True)

    ax1.set_title(r'{\bf Dome C} // $A=%.1e$, $Q=%.1e$ // $L=%i$'%(A,Q,L), fontsize=FS)

    ### Temperature
    #ax2.plot(f_T(zT), zT, 'g', lw=2)
    ax2.plot(T, z_sf, 'k', lw=2)
    ax2.set_xlabel(r'$T$ ($\SI{}{\degreeCelsius}$)')
    plt.setp(ax2.get_yticklabels(), visible=False)

    ### Gamma
    ax3.semilogx(gam, z_sf, 'k', lw=2)
    ax3.set_xlabel(r'$\Gamma_0$')
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax3.set_xticks([1e-16,1e-14,1e-12])

    ### ODFs

    import scipy.special as sp
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import cmasher as cmr
    import cartopy.crs as ccrs

    inclination, rot = 45, +45 # view angle
    prj, geo = ccrs.Orthographic(rot, 90-inclination), ccrs.Geodetic()
    axODF1 = fig.add_subplot(gs[0, 1], projection=prj)
    axODF2 = fig.add_subplot(gs[1, 1], projection=prj)
    axODF3 = fig.add_subplot(gs[2, 1], projection=prj)
    axODF1.set_global() 
    axODF2.set_global() 
    axODF3.set_global() 

    # Plot ODFs
    nn_list = [int(Nt*2/10), int(Nt*2/3), -1]
    ax_list = [axODF1, axODF2, axODF3]
    cax = '#ff7f00'

    for ii, nn in enumerate(nn_list):

        ax1.plot([0,1],[z_sf[nn],]*2, ':', c='k', lw=2)
        plot_ODF(nlm[nn,:], lm, ax=ax_list[ii], cblabel='$\psi/N$ at z=%im'%(z_sf[nn]))
        ax_list[ii].plot([0],[90], marker=r'$z$', ms=7, c=cax, transform=geo) # z axis
        ax_list[ii].plot([90],[0], marker=r'$y$', ms=7, c=cax, transform=geo) # y axis
        ax_list[ii].plot([0],[0],  marker=r'$x$', ms=7, c=cax, transform=geo) # x axis
        
    #####

    fname = 'calibrate-ddrx-L%i-%s%s.png'%(L, PROFILE, '--with-CDRX' if WITH_CDRX else '')
    print('saving %s'%(fname))
    plt.savefig(fname, dpi=200)

