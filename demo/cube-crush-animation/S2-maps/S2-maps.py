#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import copy, os, sys, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp
from scipy.spatial.transform import Rotation as R
sys.path.insert(0, '../..')
os.system('mkdir -p ./frames')

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

### Run options

MAKE_FRAME_Eij = 1
MAKE_FRAME_vi  = 1

MAKE_GIFS = 1

transparent = True

### Mode of deformation

mod = dict(type='ps', axis=2, T=1, r=0) # uniaxial compression along z
strain_target = -0.98
Nt = 100

### CPO dynamics

L = 18 # Spectral truncation
lm, nlm_len = sf.init(L)

nlm_iso = np.zeros((nlm_len), dtype=np.complex128)
nlm_iso[0] = np.sqrt(4*np.pi)

### Enhancement factors

(Eij_grain, alpha, n_grain) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)

### Phase velocities

# Physical parameters
rho = sfconst.ice['density'] # density of ice
Cij = sfconst.ice['elastic']['Bennett1968'] # Bennett (1968) parameters
lame_grain = sf.Cij_to_Lame_tranisotropic(Cij) 
alpha = 0.5 # Voigt--Reuss weight, where 0.5 = Hill average

vi_iso = sf.Vi_elastic_tranisotropic(nlm_iso, alpha, lame_grain, rho, 0,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
vP_iso, vS1_iso, vS2_iso = vi_iso[2,0], vi_iso[0,0], vi_iso[1,0]

def vi_rel(vi, vi_iso):
    return 100*(np.divide(vi,vi_iso) -1)

### Debug 

nlm_sm = np.zeros((nlm_len), dtype=np.complex128)
nlm_sm[:sf.L4len] = sf.nlm_ideal([0,0,1], 0, 4) 

vi_sm_vert = sf.Vi_elastic_tranisotropic(nlm_sm, alpha, lame_grain, rho, 0,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
vi_sm_hori = sf.Vi_elastic_tranisotropic(nlm_sm, alpha, lame_grain, rho, 90,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]

print('vi_sm_vert/vi_iso = ', vi_rel(vi_sm_vert,vi_iso))
print('vi_sm_hori/vi_iso = ', vi_rel(vi_sm_hori,vi_iso))

### Plot

geo, prj = sfplt.getprojection(rotation=50, inclination=50)

def setup_fig():

    scale = 2.6
    fig = plt.figure(figsize=(4/2*1.2*scale,0.90*scale))
    gs = fig.add_gridspec(1,4)
    al = 0.025
    ar = 0.02
    gs.update(left=al, right=1-ar, top=0.88, bottom=0.25, wspace=0.4, hspace=0.4)

    ax1 = fig.add_subplot(gs[0,0], projection=prj); ax1.set_global(); 
    ax2 = fig.add_subplot(gs[0,1], projection=prj); ax2.set_global(); 
    ax3 = fig.add_subplot(gs[0,2], projection=prj); ax3.set_global(); 
    ax4 = fig.add_subplot(gs[0,3], projection=prj); ax4.set_global(); 
    axlist = [ax1,ax2,ax3,ax4]

    fraction = 0.07
    aspect = 10

    return axlist, fig, fraction, aspect

def mkframe_Eij(nlm, Err,Ert,Erp, nn):
    
    axlist, fig, fraction, aspect = setup_fig()    
    ax1,ax2,ax3,ax4 = axlist

    # Plot ODF
    lvls = np.linspace(0.0,0.8,9)
    sfplt.plotODF(nlm, lm, ax1, lvlset=[lvls, lambda x,p:'%.1f'%x], cbaspect=aspect, cbfraction=fraction)

    # Plot Err
    cmap = 'RdBu'
    lvls = [0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 9]
    tickintvl = 2
    tick_labels = [str(v) if v <=1 else str(int(v)) for v in lvls]
    norm = colors.TwoSlopeNorm(vmin=lvls[0], vcenter=1, vmax=lvls[-1])
    plot_field(Err, ax2, lvls, cmap=cmap, cblbl=r'$E_{rr}$', tickintvl=tickintvl, norm=norm, aspect=aspect, fraction=fraction, tick_labels=tick_labels) # , norm=colors.LogNorm()

    # Plot Ert
    plot_field(Ert, ax3, lvls, cmap=cmap, cblbl=r'$E_{r\theta}$', tickintvl=tickintvl, norm=norm, aspect=aspect, fraction=fraction, tick_labels=tick_labels)

    # Plot Erp
    plot_field(Erp, ax4, lvls, cmap=cmap, cblbl=r'$E_{r\phi}$', tickintvl=tickintvl, norm=norm, aspect=aspect, fraction=fraction, tick_labels=tick_labels)

    # Axis labels
    for axi in axlist: sfplt.plotcoordaxes(axi, geo, axislabels='vuxi', color='k')
    
    # Title
    ttl = ax2.set_title(r'Uniaxial compression, $\epsilon_{zz} = %+.2f$'%(strain[nn]), fontsize=FS+2, pad=13)
    ttl.set_position([1.15, 1.07])

    # Save fig

    fout = 'frames/S2-Eij-%03i.png'%(nn)
    print('Saving %s'%(fout))
    plt.savefig(fout, transparent=transparent, dpi=200)
    plt.close()
   
def mkframe_vi(nlm, vP, vS1, vS2, nn):

    axlist, fig, fraction, aspect = setup_fig()    
    ax1,ax2,ax3,ax4 = axlist

    # Plot ODF
    lvls = np.linspace(0.0,0.8,9)
    sfplt.plotODF(nlm, lm, ax1, lvlset=[lvls, lambda x,p:'%.1f'%x], cbaspect=aspect, cbfraction=fraction)

    # Plot vP
    lvls = np.linspace(-6,6,9)
    tickintvl = 2
    plot_field(vP, ax2, lvls, cmap='RdBu', cblbl=r'$V_{P}/V^{\mathrm{iso}}_{P}-1$  (\%)', tickintvl=tickintvl, aspect=aspect, fraction=fraction) 

    # Plot vS1
    plot_field(vS1, ax3, lvls, cmap='RdBu', cblbl=r'$V_{\mathrm{S}1}/V^{\mathrm{iso}}_{\mathrm{S}1}-1$  (\%)', tickintvl=tickintvl, aspect=aspect, fraction=fraction)

    # Plot vS2
    plot_field(vS2, ax4, lvls, cmap='RdBu', cblbl=r'$V_{\mathrm{S}2}/V^{\mathrm{iso}}_{\mathrm{S}2}-1$  (\%)', tickintvl=tickintvl, aspect=aspect, fraction=fraction)

    # Axis labels
    for axi in axlist: sfplt.plotcoordaxes(axi, geo, axislabels='vuxi', color='k')

    # Title
    ttl = ax2.set_title(r'Uniaxial compression, $\epsilon_{zz} = %+.2f$'%(strain[nn]), fontsize=FS+2, pad=13)
    ttl.set_position([1.15, 1.07])

    # Save fig

    fout = 'frames/S2-vi-%03i.png'%(nn)
    print('Saving %s'%(fout))
    plt.savefig(fout, transparent=transparent, dpi=200)
    plt.close()


def plot_field(F, ax, lvls, cmap='Greys', cblbl=r'$E_{??}$', titlestr='', tickintvl=1, norm=None, tick_labels=None, aspect=9, fraction=0.6):

    hdistr = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), extend='both', norm=norm, cmap=cmap, levels=lvls, nchunk=5)
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    if not np.isscalar(lvls): 
        cb1 = plt.colorbar(hdistr, ax=ax, fraction=fraction, aspect=aspect,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])   
    else:
        cb1 = plt.colorbar(hdistr, ax=ax, fraction=fraction, aspect=aspect,  orientation='horizontal', pad=0.1)   
    cb1.set_label(cblbl)
    if tick_labels is not None: cb1.set_ticklabels(tick_labels[::tickintvl])
    ax.set_title(titlestr, fontsize=FS, pad=10)

#----- MAIN -----

if MAKE_FRAME_Eij or MAKE_FRAME_vi:

    ### Determine fabric evolution
    nlm, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=+1, nu=1) # latrot only
    strain = np.array([sf.F_to_strain(F[tt,:,:])[2,2] for tt in range(Nt+1)])

    ### S2 resolution
    latres = 40
    #latres = 25
    lat = np.deg2rad(np.linspace(-90, 90, latres))
    lon = np.deg2rad(np.linspace(0, 360, 2*latres))
    colat = np.deg2rad(90) - lat

    ### Init vars
    
    # r = radial vector, t = theta vector, p = phi vector
    dims = (Nt+1, len(lat),len(lon))
    Err = np.ones(dims) # r-r longitidinal
    Ert = np.ones(dims) # r-t shear
    Erp = np.ones(dims) # r-p shear
    vP, vS1, vS2  = np.zeros(dims), np.zeros(dims), np.zeros(dims) 

    ### Determine spherical basis vectors
    lon2, colat2 = np.meshgrid(lon, colat)
    vr, vt, vp = sfdsc.sphericalbasisvectors(colat2, lon2)

    rr = np.einsum('ikl,jkl->ijkl', vr, vr)
    rt = np.einsum('ikl,jkl->ijkl', vr, vt)
    rp = np.einsum('ikl,jkl->ijkl', vr, vp)

    id4 = np.zeros((3,3,latres,2*latres)) # repeated identity matrix
    id4[0,0,:,:] = id4[1,1,:,:] = id4[2,2,:,:] = 1
    tau_rr = id4-3*rr
    tau_rt = rt + np.einsum('ijkl->jikl',rt)
    tau_rp = rp + np.einsum('ijkl->jikl',rp)

    for nn in np.arange(0,Nt+1):

        ### Determine Eij and Vi for a given deformation (given time)

        for ii, theta in enumerate(colat):
            for jj, phi in enumerate(lon):

                r,t,p = vr[:,ii,jj], vt[:,ii,jj], vp[:,ii,jj]

                args = (Eij_grain,alpha,n_grain)
                Err[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], r,r,tau_rr[:,:,ii,jj], *args)
                Ert[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], r,t,tau_rt[:,:,ii,jj], *args)
                Erp[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], r,p,tau_rp[:,:,ii,jj], *args)

                vi = sf.Vi_elastic_tranisotropic(nlm[nn,:], alpha, lame_grain, rho, theta,phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
                vP[nn,ii,jj], vS1[nn,ii,jj], vS2[nn,ii,jj] = vi[2,0], vi[0,0], vi[1,0]
                
        ### Plot
        
        if MAKE_FRAME_Eij: mkframe_Eij(nlm[nn,:], Err[nn,:,:], Ert[nn,:,:], Erp[nn,:,:], nn)
        
        vP_rel  = vi_rel(vP[nn,:,:],  vP_iso)
        vS1_rel = vi_rel(vS1[nn,:,:], vS1_iso)
        vS2_rel = vi_rel(vS2[nn,:,:], vS2_iso)
        if MAKE_FRAME_vi:  mkframe_vi( nlm[nn,:], vP_rel, vS1_rel, vS2_rel, nn)


if MAKE_GIFS:

    for ii in np.arange(Nt,Nt+40):
        os.system('cp frames/S2-Eij-%03d.png frames/S2-Eij-%03d.png'%(Nt-1, ii))
        os.system('cp frames/S2-vi-%03d.png  frames/S2-vi-%03d.png'%(Nt-1, ii))
        
    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/S2-Eij-%03d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p S2-Eij.avi')
    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/S2-vi-%03d.png  -vcodec libx264 -crf 20  -pix_fmt yuv420p S2-vi.avi')

    os.system('rm S2-Eij.gif S2-vi.gif')

    os.system('ffmpeg -i S2-Eij.avi -vf "fps=20,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 S2-Eij.gif')
    os.system('ffmpeg -i S2-vi.avi  -vf "fps=20,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 S2-vi.gif')

       
