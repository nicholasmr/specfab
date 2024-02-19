#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import copy, os, sys, code # code.interact(local=locals())
os.system('mkdir -p ./frames')

import numpy as np
import scipy.special as sp

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc, colors
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from matplotlib.ticker import LogFormatter 

#---------------
# Run options
#---------------

MAKE_FRAME_Eij   = 0
MAKE_FRAME_vi    = 0 # elastic
MAKE_FRAME_vi_EM = 0

MAKE_GIFS = 1

transparent = True

#---------------
# Settings
#---------------

### Mode of deformation

mod = dict(type='ps', axis=2, T=1, r=1) # uniaxial compression along z
strain_target = -0.98
Nt = 100

### CPO dynamics

L = 18 # Spectral truncation
lm, nlm_len = sf.init(L)

nlm_iso = np.zeros((nlm_len), dtype=np.complex128)
nlm_iso[0] = np.sqrt(4*np.pi)

### Enhancement factors

(Eij_grain, alpha, n_grain) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)

### Elastic phase velocities

# Physical parameters
rho = sfconst.ice['density'] # density of ice
Cij = sfconst.ice['elastic']['Bennett1968'] # Bennett (1968) parameters
lame_grain = sf.Cij_to_Lame_tranisotropic(Cij) 
alpha = 0.5 # Voigt--Reuss weight, where 0.5 = Hill average

vi_iso = sf.Vi_elastic_tranisotropic(nlm_iso, alpha, lame_grain, rho, 0,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
vP_iso, vS1_iso, vS2_iso = vi_iso[2,0], vi_iso[0,0], vi_iso[1,0]

def vi_rel(vi, vi_iso):
    return 100*(np.divide(vi,vi_iso) -1)

### EM phase velocities

epsa = 3.17 - 0.034 * 1
epsc = 3.17
mu = 1

vi_EM_iso = sf.Vi_electromagnetic_tranisotropic(nlm_iso, epsc,epsa,mu, 0,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:]
vS1_EM_iso, vS2_EM_iso = vi_EM_iso[0,0], vi_EM_iso[1,0]
print('vSi EM:',  vS1_EM_iso*1e-6, vS2_EM_iso*1e-6 )

def vi_EM_rel(vi, vi_iso):
    return 100*(np.divide(vi,vi_iso) -1)

### Debug 

nlm_sm = np.zeros((nlm_len), dtype=np.complex128)
nlm_sm[:sf.L4len] = sf.nlm_ideal([0,0,1], 0, 4) 

vi_sm_vert = sf.Vi_elastic_tranisotropic(nlm_sm, alpha, lame_grain, rho, 0,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
vi_sm_hori = sf.Vi_elastic_tranisotropic(nlm_sm, alpha, lame_grain, rho, 90,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]

print('vi_sm_vert/vi_iso = ', vi_rel(vi_sm_vert,vi_iso))
print('vi_sm_hori/vi_iso = ', vi_rel(vi_sm_hori,vi_iso))

#---------------
# Plotting routines
#---------------

geo, prj = sfplt.getprojection(rotation=50, inclination=50)

def setup_fig():

    scale = 2.6
    fig = plt.figure(figsize=(4/2*1.2*scale,0.90*scale))
    gs = fig.add_gridspec(1,4)
    gs.update(left=0.025, right=1-0.02, top=0.88, bottom=0.25, wspace=0.4, hspace=0.4)

    axlist = [fig.add_subplot(gs[0,ii], projection=prj) for ii in range(4)]
    for axi in axlist: axi.set_global()

    fraction = 0.07
    aspect = 10

    return axlist, fig, fraction, aspect


def setup_fig_narrow():

    scale = 2.6
    fig = plt.figure(figsize=(3/2*1.1*scale,0.90*scale))
    gs = fig.add_gridspec(1,3)
    gs.update(left=0.025, right=1-0.02, top=0.88, bottom=0.25, wspace=0.4, hspace=0.4)

    axlist = [fig.add_subplot(gs[0,ii], projection=prj) for ii in range(3)]
    for axi in axlist: axi.set_global()

    fraction = 0.07
    aspect = 10

    return axlist, fig, fraction, aspect
    

def mkframe_Eij(nlm, Err,Ert,Erp, mlon,mlat, nn):
    
    axlist, fig, fraction, aspect = setup_fig()    
    ax1,ax2,ax3,ax4 = axlist

    sfplt.plotODF(nlm, lm, ax1, lvlset=[np.linspace(0.0,0.8,9), lambda x,p:'%.1f'%x], cbaspect=aspect, cbfraction=fraction)

    lvls = [0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 9]
    tick_labels = [str(v) if v <=1 else str(int(v)) for v in lvls]
    norm = colors.TwoSlopeNorm(vmin=lvls[0], vcenter=1, vmax=lvls[-1])
    kwargs = dict(cmap='RdBu', tickintvl=2, norm=norm, aspect=aspect, fraction=fraction, tick_labels=tick_labels)

    plot_field(Err, mlon, mlat, ax2, lvls, cblbl=r'$E_{rr}$', **kwargs)
    plot_field(Ert, mlon, mlat, ax3, lvls, cblbl=r'$E_{r\theta}$', **kwargs)
    plot_field(Erp, mlon, mlat, ax4, lvls, cblbl=r'$E_{r\phi}$', **kwargs)

    for axi in axlist: sfplt.plotcoordaxes(axi, geo, axislabels='vuxi', color='k')
    
    ttl = ax2.set_title(r'Uniaxial compression, $\epsilon_{zz} = %+.2f$'%(strain[nn]), fontsize=FS+2, pad=13)
    ttl.set_position([1.15, 1.07])

    fout = 'frames/S2-Eij-%03i.png'%(nn)
    print('Saving %s'%(fout))
    plt.savefig(fout, transparent=transparent, dpi=200)
    plt.close()
   
   
def mkframe_vi(nlm, vP,vS1,vS2, mlon,mlat, nn):

    axlist, fig, fraction, aspect = setup_fig()    
    ax1,ax2,ax3,ax4 = axlist

    sfplt.plotODF(nlm, lm, ax1, lvlset=[np.linspace(0.0,0.8,9), lambda x,p:'%.1f'%x], cbaspect=aspect, cbfraction=fraction)

    kwargs = dict(cmap='RdBu', tickintvl=2, aspect=aspect, fraction=fraction)
    lvls = np.linspace(-6,6,9)

    plot_field(vP,  mlon, mlat, ax2, lvls, cblbl=r'$V_{P}/V^{\mathrm{iso}}_{P}-1$  (\%)', **kwargs)
    plot_field(vS1, mlon, mlat, ax3, lvls, cblbl=r'$V_{\mathrm{S}1}/V^{\mathrm{iso}}_{\mathrm{S}1}-1$  (\%)', **kwargs)
    plot_field(vS2, mlon, mlat, ax4, lvls, cblbl=r'$V_{\mathrm{S}2}/V^{\mathrm{iso}}_{\mathrm{S}2}-1$  (\%)', **kwargs)

    for axi in axlist: sfplt.plotcoordaxes(axi, geo, axislabels='vuxi', color='k')

    ttl = ax2.set_title(r'Uniaxial compression, $\epsilon_{zz} = %+.2f$'%(strain[nn]), fontsize=FS+2, pad=13)
    ttl.set_position([1.15, 1.07])

    fout = 'frames/S2-vi-%03i.png'%(nn)
    print('Saving %s'%(fout))
    plt.savefig(fout, transparent=transparent, dpi=200)
    plt.close()


def mkframe_vi_EM(nlm, vS1,vS2,vSrel, mlon,mlat, nn):

    axlist, fig, fraction, aspect = setup_fig_narrow()
    ax1,ax2,ax3 = axlist
#    axlist, fig, fraction, aspect = setup_fig()    
#    ax1,ax2,ax3,ax4 = axlist

    sfplt.plotODF(nlm, lm, ax1, lvlset=[np.linspace(0.0,0.8,9), lambda x,p:'%.1f'%x], cbaspect=aspect, cbfraction=fraction)

    kwargs = dict(cmap='RdBu', tickintvl=4, aspect=aspect, fraction=fraction)
    lvls = np.linspace(-0.2,0.2,9)

    plot_field(vS1, mlon, mlat, ax2, lvls, cblbl=r'$V_{\mathrm{S}1}/V^{\mathrm{iso}}_{\mathrm{S}1}-1$  (\%)', **kwargs)
    plot_field(vS2, mlon, mlat, ax3, lvls, cblbl=r'$V_{\mathrm{S}2}/V^{\mathrm{iso}}_{\mathrm{S}2}-1$  (\%)', **kwargs)
    
    #lvls = np.linspace(-2,2,9)
    #plot_field(vSrel, mlon, mlat, ax4, lvls, cblbl=r'$V_{\mathrm{S}1}/V_{\mathrm{S}1}-1$  (\%)', **kwargs)

    for axi in axlist: sfplt.plotcoordaxes(axi, geo, axislabels='vuxi', color='k')

    ttl = ax2.set_title(r'Confined compression, $\epsilon_{zz} = %+.2f$'%(strain[nn]), fontsize=FS+2, pad=13)
    ttl.set_position([0.5, 1.07])

    fout = 'frames/S2-vi-EM-%03i.png'%(nn)
    print('Saving %s'%(fout))
    plt.savefig(fout, transparent=transparent, dpi=200)
    plt.close()


def plot_field(F, mlon, mlat, ax, lvls, cmap='Greys', cblbl=r'$E_{??}$', titlestr='', tickintvl=1, norm=None, tick_labels=None, aspect=9, fraction=0.6):

    kwargs_cf = dict(extend='both', norm=norm, cmap=cmap, levels=lvls, nchunk=5)
    kwargs_cb = dict(ticks=lvls[::tickintvl], fraction=fraction, aspect=aspect, orientation='horizontal', pad=0.1)    
    hcf, hcb = sfplt.plotS2field(ax, F, mlon, mlat, kwargs_cf=kwargs_cf, showcb=True, kwargs_cb=kwargs_cb)
    hcb.set_label(cblbl)
    if tick_labels is not None: hcb.set_ticklabels(tick_labels[::tickintvl])
    ax.set_title(titlestr, fontsize=FS, pad=10)
    

#---------------
# Main
#---------------

if MAKE_FRAME_Eij or MAKE_FRAME_vi or MAKE_FRAME_vi_EM:

    ### Determine fabric evolution
    
    nlm, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=+1, nu=1) # latrot only
    strain = np.array([sf.F_to_strain(F[tt,:,:])[2,2] for tt in range(Nt+1)])

    ### Init vars    

    # r = radial vector, t = theta vector, p = phi vector
    (vr,vt,vp, rr,rt,rp, tau_rr,tau_rt,tau_rp, vlon,vlat, mlon,mlat) = sfdsc.idealstressstate()

    dims = (Nt+1, len(vlat),len(vlon))
    Err = np.ones(dims) # r-r longitudinal
    Ert = np.ones(dims) # r-t shear
    Erp = np.ones(dims) # r-p shear
    vP, vS1, vS2 = np.zeros(dims), np.zeros(dims), np.zeros(dims) 
    vS1_EM, vS2_EM = np.zeros(dims), np.zeros(dims)

    ### Determine Eij and Vi for a given deformation (given time)

    for nn in np.arange(0,Nt+1):

        for ii, theta in enumerate(vlat):
            for jj, phi in enumerate(vlon):

                colat = sfdsc.lat2colat(theta)

                r,t,p = vr[:,ii,jj], vt[:,ii,jj], vp[:,ii,jj]

                if MAKE_FRAME_Eij: 
                    args = (Eij_grain,alpha,n_grain)
                    Err[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], r,r,tau_rr[:,:,ii,jj], *args)
                    Ert[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], r,t,tau_rt[:,:,ii,jj], *args)
                    Erp[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], r,p,tau_rp[:,:,ii,jj], *args)

                if MAKE_FRAME_vi: 
                    vi = sf.Vi_elastic_tranisotropic(nlm[nn,:], alpha, lame_grain, rho, colat,phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
                    vP[nn,ii,jj], vS1[nn,ii,jj], vS2[nn,ii,jj] = vi[2,0], vi[0,0], vi[1,0]
                
                if MAKE_FRAME_vi_EM: 
                    vS1_EM[nn,ii,jj], vS2_EM[nn,ii,jj] = sf.Vi_electromagnetic_tranisotropic(nlm[nn,:], epsc,epsa,mu, colat,phi) 
               
        ### Plot
        
        if MAKE_FRAME_Eij: 
            mkframe_Eij(nlm[nn,:], Err[nn,:,:], Ert[nn,:,:], Erp[nn,:,:], mlon, mlat, nn)
        
        if MAKE_FRAME_vi: 
            vP_rel  = vi_rel(vP[nn,:,:],  vP_iso)
            vS1_rel = vi_rel(vS1[nn,:,:], vS1_iso)
            vS2_rel = vi_rel(vS2[nn,:,:], vS2_iso)
            mkframe_vi(nlm[nn,:], vP_rel, vS1_rel, vS2_rel, mlon, mlat, nn)

        if MAKE_FRAME_vi_EM: 
            vS1_EM_rel = vi_EM_rel(vS1_EM[nn,:,:], vS1_EM_iso)
            vS2_EM_rel = vi_EM_rel(vS2_EM[nn,:,:], vS2_EM_iso)
            vS_EM_rel = np.divide(vS1_EM[nn,:,:],vS2_EM[nn,:,:])
            mkframe_vi_EM(nlm[nn,:], vS1_EM_rel,vS2_EM_rel,vS_EM_rel, mlon, mlat, nn)


if MAKE_GIFS:

    for ii in np.arange(Nt,Nt+40):
        os.system('cp frames/S2-Eij-%03d.png frames/S2-Eij-%03d.png'%(Nt-1, ii))
        os.system('cp frames/S2-vi-%03d.png  frames/S2-vi-%03d.png'%(Nt-1, ii))
        os.system('cp frames/S2-vi-EM-%03d.png  frames/S2-vi-EM-%03d.png'%(Nt-1, ii))
        
    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/S2-Eij-%03d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p S2-Eij.avi')
    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/S2-vi-%03d.png  -vcodec libx264 -crf 20  -pix_fmt yuv420p S2-vi.avi')
    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/S2-vi-EM-%03d.png  -vcodec libx264 -crf 20  -pix_fmt yuv420p S2-vi-EM.avi')

    os.system('rm S2-Eij.gif S2-vi.gif S2-vi-EM.gif')

    os.system('ffmpeg -i S2-Eij.avi   -vf "fps=20,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 S2-Eij.gif')
    os.system('ffmpeg -i S2-vi.avi    -vf "fps=20,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 S2-vi.gif')
    os.system('ffmpeg -i S2-vi-EM.avi -vf "fps=20,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 S2-vi-EM.gif')

       
