# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2022-2023

# Calibrate DDRX activation function energy (Q) and prefactor (A) by reproducing Dome C (EDC) fabric profile.

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from scipy.interpolate import interp1d
import pickle
from progress.bar import Bar

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

yr2s = 31556926
s2yr = 3.16887646e-8
c2k  = lambda degc: degc+273.15 # deg. C to deg. K

#---------------------
# Numerics
#---------------------

L_list = [4,6,8,20] # Spectral truncation used for calibration
L_list = [8,]

Nt = 200 # Number of integration steps

WITH_CDRX = 0

for L in L_list:

    #---------------------
    # Gamma_0 
    #---------------------

    Q = 1e5 # Activation energy

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
        file = open("EDC.p",'rb')
        (fab,temp,aux) = pickle.load(file)
        n20_0 = 0.14 # initial fabric state of parcel (only l,m=2,0 component is assumed nonzero)
        strain_zz_stop=-0.96 # stop integration at this vertical parcel strain. Fabric is only measured until the depth corresponding to this vertical strain.

    # Geometry
    z_bed = -aux['H'] # depth of bed
    z0 = abs(fab['z'][0]) 
    z1 = aux['H']
    H = z1-z0 # ice column height as seen by parcel model
    print('z_bed, H = %f, %f'%(z_bed, H))

    # Measured fabric
    lami_meas, z_meas = fab['eig'], fab['z']

    # Temperature interpolant
    T, zT = temp['T'], temp['z']
    f_T = interp1d(zT,T) # C. Ritz (pers. comm.)

    # Accumulation rate
    #b = aux['b'] # present day
    b = 0.0153 # average over core (D. Lilien)
    print('b = %3e m/yr'%(b))

    #---------------------
    # Ice parcel model
    #---------------------
        
    # Assumes depth-constant (unconfined) vertical compression: the classical Nye model of a dome.
    b /= yr2s
    t_e = 1/(b/H) # e-folding time for uniaxial compression (strainrate_zz = MeanAccum/H)
    r = 0 # uniaxial compression
    ax = 2 # z axis
        
    # Constants
    tend    = sf.pureshear_strainii_to_t(strain_zz_stop, t_e)
    timevec = np.linspace(0,tend,Nt)
    dt      = timevec[1]-timevec[0] # given t_e, if we want Nt time steps, this is the step size needed
    strain_zz = np.array([sf.F_to_strain(sf.pureshear_F(ax,r,t_e,t))[-1,-1] for ii, t in enumerate(timevec)]) # for the constant time-step size, these are the vertical parcel strains as a function of time
    z_sf = H*(strain_zz+1) - H # Vertical parcel strains correspond to these depths
    z_sf += -z0 # effective offset of parcel model

    # Initialize structure
    lm, nlm_len = sf.init(L)
    nlm = np.zeros((Nt, nlm_len), dtype=np.complex64) # array of expansion coefficients
    nlm[0,0] = 1/np.sqrt(4*np.pi) # normalized distribution
    nlm[0,3] = n20_0 # initial single max strength (manually fitted to match eigenvalue profile)
    lami_sf = np.zeros((Nt, 3)) # array of modelled eigen values
    lami_sf[0,:] = np.diag(sf.a2(nlm[0,:])) # a2 already diagional for this mode of deformation, so simply extract the diagonal values
    #
    gam = np.zeros((Nt)) # DDRX rate factor magnitude
    lam = np.zeros((Nt)) # CDRX rate factor magnitude
        
    # Euler integration of fabric model
    D, W = sf.ugrad_to_D_and_W(sf.pureshear_ugrad(ax,r,t_e)) # strain-rate and spin tensors for mode of deformation
    T = f_T(z_sf) # Temperature at modelled parcel depths 
    
    def f_lam(D,T):
        m_lam, b_lam = 1.26e-3, 0.206 # Richards et al. (2021), Lilien et al. (2022) values
        eps_E = np.sqrt(0.5*np.einsum('ij,ji',D,D)) # sqrt(0.5 * D:D)
        return eps_E * (m_lam*T + b_lam)

    lam[0] = f_lam(D,T[0])
    gam[0] = sf.Gamma0(D, c2k(T[0]), A, Q)
    
    with Bar('dt=%.3fyr, Nt=%i :: L=%i (nlm_len=%i) ::'%(dt*s2yr,Nt,L,nlm_len), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds') as bar:
        for nn in np.arange(1,Nt):
        
            nlm_prev = nlm[nn-1,:]
            
            # Lattice rotation operator (nlm_len x nlm_len matrix)
            M_LROT = sf.M_LROT(nlm_prev, D, W, 1, 0) 

            # DDRX operator (nlm_len x nlm_len matrix)
            S = D.copy() # S (dev. stress tensor) remains coaxial with D (strain-rate tensor) for this mode of deformation, so we need not calculate S from the bulk flow law (which in turn requires calculating the enhancement factors from the modelled fabric)
            gam[nn] = sf.Gamma0(D, c2k(T[nn]), A, Q) # DDRX decay rate magnitude
            M_DDRX = gam[nn] * sf.M_DDRX(nlm_prev, S) # DDRX operator

            # CDRX operator (nlm_len x nlm_len matrix)
            if WITH_CDRX:
                lam[nn] = f_lam(D,T[nn]) 
                M_CDRX = lam[nn]*sf.M_CDRX(nlm_prev) 
            else:
                M_CDRX = 0*M_DDRX

            # Regularization operator (nlm_len x nlm_len matrix)
            M_REG = sf.M_REG(nlm_prev, D) 

            # Total fabric evolution 
            M = M_LROT + M_DDRX + M_CDRX + M_REG # net operator
            nlm[nn,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # Forward Euler step
            nlm[nn,:] = sf.apply_bounds(nlm[nn,:]) # Apply spectral bounds, if needed 
            
            lami_sf[nn,:] = np.diag(sf.a2(nlm[nn,:])) # a2 already diagional for this mode of deformation, so simply extract the diagonal values
            
            bar.next() # update progress bar
            
    # Estimate numerical accuracy 
    diff_abs = np.abs(100*(nlm[:,0]-nlm[0,0])/nlm[0,0])
    print(r'Numerical error estimate: max(n00(t) - n00(t=0)) [%] = ', np.amax(diff_abs))

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

    geo, prj = sfplt.getprojection(rotation=45, inclination=45)
    
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
        sfplt.plotODF(nlm[nn,:], lm, ax_list[ii], lvlset='zero-up', cblabel='$n/N$ at z=%im'%(z_sf[nn]))
        sfplt.plotcoordaxes(ax_list[ii], geo, axislabels='vuxi')
        
    #####

    fname = 'calibrate-ddrx-L%i-%s%s.png'%(L, PROFILE, '--with-CDRX' if WITH_CDRX else '')
    print('saving %s'%(fname))
    plt.savefig(fname, dpi=150)

