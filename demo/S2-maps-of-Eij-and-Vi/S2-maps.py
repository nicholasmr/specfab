#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import copy, os, sys, code # code.interact(local=locals())
import numpy as np
sys.path.insert(0, '..')
os.system('mkdir -p ./frames')

from specfabpy import specfabpy as sf
from header import *
from sfconstants import *

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc, colors
import matplotlib.gridspec as gridspec
import scipy.special as sp
from matplotlib.ticker import LogFormatter 

MAKE_FRAME_Eij = 1
MAKE_FRAME_vi  = 1

MAKE_GIFS = 1

transparent = False

### CPO dynamics

f = 2
Nt = f*40 # Number of time steps
dt = 1/f*0.0782404601085629 # Time-step size (gives a vertical strain of -0.98 for experiment "uc_zz")
ugrad = np.diag([0.5, 0.5, -1.0]) # Uniaxial compression along z
eps = (ugrad+np.transpose(ugrad))/2 # Symmetric part (strain-rate)
omg = (ugrad-np.transpose(ugrad))/2 # Anti-symmetric part (spin)
te = 1/eps[-1,-1]
strainzz = lambda t: np.exp(t/te)-1
print(r'Integrating until strain_zz = %.3f'%(strainzz(Nt*dt)))

#L = 12 # Spectral truncation
L = 18
lm, nlm_len = sf.init(L) # nlm_len is the number of fabric expansion coefficients (degrees of freedom).
nlm = np.zeros((Nt,nlm_len), dtype=np.complex128)
nlm[:,0] = 1/np.sqrt(4*np.pi) # Normalized such that N(t=0) = 1

### Enhancement factors

(Eij_grain, alpha, n_grain) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)

latres = 50
lat = np.deg2rad(np.linspace(-90,90,latres))
lon = np.deg2rad(np.linspace(0,360,2*latres))
colat = np.deg2rad(90) - lat

# r = radial vector, t = theta vector, p = phi vector
dims = (Nt, len(lat),len(lon))
Err = np.ones(dims) # r-r longitidinal
Ert = np.ones(dims) # r-t shear
Erp = np.ones(dims) # r-p shear

### Phase velocities

# Physical parameters
rho = 917 # density of ice
C11,C33,C55,C12,C13 = 14.060e9, 15.240e9, 3.060e9, 7.150e9, 5.880e9 # Bennett (1968) parameters
lam,mu,Elam,Emu,Egam = sf.Cij_to_Lame_tranisotropic(C11,C33,C55,C12,C13) 
alpha = 0.5 # Voigt--Reuss weight, where 0.5 = Hill average

vi_iso = sf.Vi_elastic_tranisotropic(nlm[0,:], alpha, lam,mu,Elam,Emu,Egam, rho, 0,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
vP_iso, vS1_iso, vS2_iso = vi_iso[2,0], vi_iso[0,0], vi_iso[1,0]
vP, vS1, vS2  = np.zeros(dims), np.zeros(dims), np.zeros(dims) 

def vi_rel(vi, vi_iso):
    return 100*(np.divide(vi,vi_iso) -1)

### debug 

m = np.array([0,0,1])
a2_sm = np.einsum('i,j',m,m)
a4_sm = np.einsum('i,j,k,l',m,m,m,m)
n2m = sf.a2_to_nlm(a2_sm)
n4m = sf.a4_to_nlm(a4_sm)
nlm_sm = np.zeros((nlm_len), dtype=np.complex128)
#nlm[:(1+5)] = n2m[:]
nlm_sm[:(1+5+9)] = n4m[:]

vi_sm_vert = sf.Vi_elastic_tranisotropic(nlm_sm, alpha, lam,mu,Elam,Emu,Egam, rho, 0,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
vi_sm_hori = sf.Vi_elastic_tranisotropic(nlm_sm, alpha, lam,mu,Elam,Emu,Egam, rho, 90,0) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]

print('vi_sm_vert/vi_iso = ', vi_rel(vi_sm_vert,vi_iso))
print('vi_sm_hori/vi_iso = ', vi_rel(vi_sm_hori,vi_iso))

### Plot

inclination = 45 # view angle
rot0 = 1 * -90 
rot = rot0 - 45 # view angle

prj = ccrs.Orthographic(rot, 90-inclination)
geo = ccrs.Geodetic()

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

def set_axis_labels(axlist):
    FSAX = FS+1
    for ax in axlist:
        ax.text(rot0-40, 85, r'$\vu{z}$', horizontalalignment='left', transform=geo, fontsize=FSAX)
        ax.text(rot0-96, -8, r'$\vu{x}$', horizontalalignment='left', transform=geo, fontsize=FSAX)
        ax.text(rot0-3, -5, r'$\vu{y}$', horizontalalignment='left', transform=geo, fontsize=FSAX)

def mkframe_Eij(nlm, Err,Ert,Erp, nn):
    
    axlist, fig, fraction, aspect = setup_fig()    
    ax1,ax2,ax3,ax4 = axlist

    # Plot ODF
    lvls = np.linspace(0.0,0.8,9)
    tickintvl = 4
    plot_ODF(nlm,lm, ax=ax1, cmap='Greys', cblabel='$n/N$ (ODF)', lvls=lvls, tickintvl=tickintvl, aspect=aspect, fraction=fraction)

    # Plot Err
    cmap = 'RdBu'
    lvls = [0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 9]
    tickintvl = 2
    tick_labels = [str(v) if v <=1 else str(int(v)) for v in lvls]
    norm = colors.TwoSlopeNorm(vmin=lvls[0], vcenter=1, vmax=lvls[-1])
    plot_field(Err, ax2, lvls, cmap=cmap, cblbl=r'$E_{rr}$', tickintvl=tickintvl, norm=norm, aspect=aspect, fraction=fraction, tick_labels=tick_labels) # , norm=colors.LogNorm()

    # Plot Ert
#    lvls = np.linspace(1, 9, 5)
#    tickintvl = 2
    plot_field(Ert, ax3, lvls, cmap=cmap, cblbl=r'$E_{r\theta}$', tickintvl=tickintvl, norm=norm, aspect=aspect, fraction=fraction, tick_labels=tick_labels)

    # Plot Erp
    plot_field(Erp, ax4, lvls, cmap=cmap, cblbl=r'$E_{r\phi}$', tickintvl=tickintvl, norm=norm, aspect=aspect, fraction=fraction, tick_labels=tick_labels)

    # Axis labels
    set_axis_labels(axlist)
    
    # Title
    ttl = ax2.set_title(r'Uniaxial compression, $\epsilon_{zz} = %+.2f$'%(strainzz(nn*dt)), fontsize=FS+2, pad=13)
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
    tickintvl = 4
    plot_ODF(nlm,lm, ax=ax1, cmap='Greys', cblabel='$n/N$ (ODF)', lvls=lvls, tickintvl=tickintvl, aspect=aspect, fraction=fraction)

    # Plot vP
    lvls = np.linspace(-6,6,9)
    tickintvl = 2
    plot_field(vP, ax2, lvls, cmap='RdBu', cblbl=r'$V_{P}/V^{\mathrm{iso}}_{P}-1$  (\%)', tickintvl=tickintvl, aspect=aspect, fraction=fraction) 

    # Plot vS1
    plot_field(vS1, ax3, lvls, cmap='PRGn', cblbl=r'$V_{\mathrm{S}1}/V^{\mathrm{iso}}_{\mathrm{S}1}-1$  (\%)', tickintvl=tickintvl, aspect=aspect, fraction=fraction)

    # Plot vS2
    plot_field(vS2, ax4, lvls, cmap='PRGn', cblbl=r'$V_{\mathrm{S}2}/V^{\mathrm{iso}}_{\mathrm{S}2}-1$  (\%)', tickintvl=tickintvl, aspect=aspect, fraction=fraction)

    # Axis labels
    set_axis_labels(axlist)

    # Title
    ttl = ax2.set_title(r'Uniaxial compression, $\epsilon_{zz} = %+.2f$'%(strainzz(nn*dt)), fontsize=FS+2, pad=13)
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


if MAKE_FRAME_Eij or MAKE_FRAME_vi:

    nn = 0
    if MAKE_FRAME_Eij: mkframe_Eij(nlm[nn,:], Err[nn,:,:], Ert[nn,:,:], Erp[nn,:,:], nn)
    if MAKE_FRAME_vi:  mkframe_vi(nlm[nn,:], vP[nn,:,:], vS1[nn,:,:], vS2[nn,:,:], nn)

    for nn in np.arange(1,Nt):

        ### Update fabric

        nlm_prev = nlm[nn-1,:]
        M_LROT = sf.M_LROT(nlm_prev, eps, omg, 1, 0) # Here we consider constant large-scale velocity gradients, **but if for any practical time-varying scenario M_ij should be calculated inside the loop below!**
        if L > 12:
            nu0, expo = 10, 2
            M_REG = M_REG_custom(nu0, expo, eps, sf) # from header.py
        else:
            M_REG  = sf.M_REG(nlm_prev, eps)
        M = M_LROT + M_REG
        nlm[nn,:] = nlm_prev + dt*np.matmul(M, nlm_prev)

        ### Save S2 maps of Eij 
        
        for ii, theta in enumerate(lat):
            for jj, phi in enumerate(lon):

                ct, st = np.cos(theta), np.sin(theta)
                cp, sp = np.cos(phi), np.sin(phi)
            
                vr = np.array([ cp*st, sp*st, ct ])
                vt = np.array([ cp*ct, sp*ct, -st ])
                vp = np.array([ -sp, cp, 0 ])
        #        print(np.dot(r,t), np.dot(r,p), np.linalg.norm(r), np.linalg.norm(t), np.linalg.norm(p))
            
                rr = np.tensordot(vr,vr, axes=0)
                rt = np.tensordot(vr,vt, axes=0)
                rp = np.tensordot(vr,vp, axes=0)

                tau_rr = 1*(np.identity(3)-3*rr) 
                tau_rt = 1*(rt + np.transpose(rt)) 
                tau_rp = 1*(rp + np.transpose(rp)) 
                
                Err[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], vr,vr,tau_rr, Eij_grain,alpha,n_grain)
                Ert[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], vr,vt,tau_rt, Eij_grain,alpha,n_grain)
                Erp[nn,ii,jj] = sf.Evw_tranisotropic(nlm[nn,:], vr,vp,tau_rp, Eij_grain,alpha,n_grain)
        
                vi = sf.Vi_elastic_tranisotropic(nlm[nn,:], alpha, lam,mu,Elam,Emu,Egam, rho, colat[ii],phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
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
    os.system('ffmpeg -i S2-Eij.avi -vf "fps=22,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 S2-Eij.gif')
    os.system('ffmpeg -i S2-vi.avi  -vf "fps=22,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 S2-vi.gif')
       
