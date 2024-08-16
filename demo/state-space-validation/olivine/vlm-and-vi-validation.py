# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import copy, sys, code # code.interact(local=locals())
import numpy as np
import pickle
import pandas as pd
from progress.bar import Bar

sys.path.insert(0, '..')
from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
from specfabpy import constants as sfconst
from specfabpy import integrator as sfint
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-1

lm, nlm_len = sf.init(4)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cmasher as cmr

geo, prj = sfplt.getprojection(rotation=45-180, inclination=50)

### Parameters

rho = sfconst.olivine['density'] 
Cij = sfconst.olivine['elastic']['Abramson1997'] 
lame_grain = sf.Cij_to_Lame_orthotropic(Cij) 
alpha = 1 # Voigt--Reuss weight

viscparams_grain = sfconst.olivine['viscoplastic']['linear'] 

### Setup

DO_DREX = 1 # D-Rex runs (to be loaded, simulated elsewhere)
DO_DDM  = 1 # Discrete Directors Method 
#DO_SDM  = 0 # Spectral Directors Method *** NOT IMPLEMENTED (NOT NEEDED FOR ANALYSIS) ***

PLOT_vlm = 1
PLOT_Vi  = 1
PLOT_Eij = 1

# Experiments to process
deformationTypes = ['simpleShear','axisymmetricCompression']
#deformationTypes = ['simpleShear',]
deformationTypes = ['axisymmetricCompression',]

#stepsToPlot = [1,15,20]
stepsToPlot = [15,] # only 15 is used for papar plots
stepsToPlot = [0,] # debug, single crystal behaviour

# r = radial vector, t = theta vector, p = phi vector
(vr,vt,vp, rr,rt,rp, tau_rr,tau_rt,tau_rp, vlon,vlat, mlon,mlat) = sfdsc.idealstressstate(lonres=80, latres=40)
mcolat = sfdsc.lat2colat(mlat)

### Functions

lvlset = [np.linspace(0,0.2,9), lambda x,p:'%.1f'%x] # [np.linspace(0.0, 0.2, 5), lambda x,p:'%.1f'%x]

def plot_vlm_eval(blm,nlm,vlm,vlm_est, ii, deformationType):
   
    size = 2.0
    fig = plt.figure(figsize=(2.7*size,0.95*size))
    gs = fig.add_gridspec(ncols=4, nrows=1)
    axb, axn, ax1, ax2 = [fig.add_subplot(gs[0,ii], projection=prj) for ii in range(4)] 

    kwargs = dict(lvlset=lvlset, cbtickintvl=4, nchunk=0)
    sfplt.plotODF(blm, lm, axb, cblabel=r'$b/N$', cmap='Reds', **kwargs)
    sfplt.plotODF(nlm, lm, axn, cblabel=r'$n/N$', cmap='Blues', **kwargs)
    sfplt.plotODF(vlm, lm, ax1, cblabel=r'$v/N$', cmap='Greens', **kwargs)
    sfplt.plotODF(vlm_est, lm, ax2, cblabel=r'$v/N$', cmap='Greens', **kwargs)

    for ax in (axb,axn,ax1,ax2):
        ax.set_global()
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', color='k')
    
    plt.tight_layout()
    fname = 'vlm-and-vi-validation-plots/%s-%i.pdf'%(deformationType, ii)
    plt.savefig(fname, transparent=True)
    os.system('pdfcrop %s %s'%(fname,fname))

def plot_vi_eval(blm, vP,vS1,vS2, lon, lat, ii, deformationType):

    size = 2.0
    fig = plt.figure(figsize=(3/4*2.7*size,0.95*size))
    gs = fig.add_gridspec(ncols=3, nrows=1)
    ax1, ax2, ax3 = [fig.add_subplot(gs[0,ii], projection=prj) for ii in range(3)] 

    lvlsP = np.linspace(7.4, 9.2, 7)
    lvlsS1 = np.linspace(4.6, 5.2, 7) # np.linspace(4.7, 5.3, 9)
    lvlsS2 = lvlsS1
    plot_vi(vP*1e-3,  lon,lat, ax=ax1, cblabel='$V_{\mathrm{P}}$ (\SI{}{\kilo\metre\per\second})',  lvls=lvlsP)
    plot_vi(vS1*1e-3, lon,lat, ax=ax2, cblabel='$V_{\mathrm{S1}}$ (\SI{}{\kilo\metre\per\second})', lvls=lvlsS1)
    plot_vi(vS2*1e-3, lon,lat, ax=ax3, cblabel='$V_{\mathrm{S2}}$ (\SI{}{\kilo\metre\per\second})', lvls=lvlsS2)

#    lvlsdS = np.linspace(0,0.6, 9)
#    plot_vi((vS1-vS2)*1e-3, lon,lat, ax=ax4, cblabel='$V_{\mathrm{S1}}-V_{\mathrm{S2}}$ (\SI{}{\kilo\metre\per\second})', lvls=lvlsdS, cunder='tab:red')
    
    for ax in (ax1,ax2,ax3):
        ax.set_global()
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', color='k')
    
    plt.tight_layout()
    fname = 'vlm-and-vi-validation-plots/%s-%i.pdf'%(deformationType, ii)
    plt.savefig(fname, transparent=True)
    os.system('pdfcrop %s %s'%(fname,fname))


def plot_vi(F, lon,lat, ax=None, cmap='Greys', cunder=None, cblabel='$v$', rot0=-40, lvls=5, tickintvl=3, nchunk=0):
    
    if cmap=='Greys': cmap = cmr.get_sub_cmap(cmap, 0.05, 1) # don't include pure white.
    if cunder is not None: cmap.set_under(cunder)
    kwargs_cf = dict(levels=lvls, extend='both', cmap=cmap, nchunk=nchunk)
    kwargs_cb = dict(fraction=0.075, aspect=9,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])    
    hcf, hcb = sfplt.plotS2field(ax, F, mlon, mlat, kwargs_cf=kwargs_cf, showcb=True, kwargs_cb=kwargs_cb)
    hcb.set_label(cblabel)
    hcb.ax.xaxis.set_ticks(lvls, minor=True)
    return hcf
    

def plot_Eij_eval(blm,nlm,vlm, Err,Ert,Erp, lon, lat, ii, deformationType):
    
    size = 2.0
    fig = plt.figure(figsize=(3/4*2.7*size,0.95*size))
    gs = fig.add_gridspec(ncols=3, nrows=1)
    ax1, ax2, ax3 = [fig.add_subplot(gs[0,ii], projection=prj) for ii in range(3)] 

#    sfplt.plotODF(nlm, lm, ax1, cblabel=r'$n/N$', cmap='Blues',  lvlset=lvlset, nchunk=None)

    lvls = [0.4, 0.6, 0.8, 1, 1+1/3, 1+2/3, 2]
#    tick_labels = [str(v) if v <=1 else str(int(v)) for v in lvls]
    tick_labels = [v for v in lvls]
    kwargs = dict(cmap='PuOr', tickintvl=3, norm=colors.TwoSlopeNorm(vmin=lvls[0], vcenter=1, vmax=lvls[-1]), aspect=10, fraction=0.07, tick_labels=tick_labels)
    plot_Eij(Err, lon, lat, ax1, lvls, cblbl=r'$E_{rr}$', **kwargs)
    plot_Eij(Ert, lon, lat, ax2, lvls, cblbl=r'$E_{r\theta}$', **kwargs)
    plot_Eij(Erp, lon, lat, ax3, lvls, cblbl=r'$E_{r\phi}$', **kwargs)

    for ax in (ax1,ax2,ax3):
        ax.set_global()
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', color='k')
        
    plt.tight_layout()
    fname = 'vlm-and-vi-validation-plots/%s-%i.pdf'%(deformationType, ii)
    plt.savefig(fname, transparent=True)
    os.system('pdfcrop %s %s'%(fname,fname))


def plot_Eij(F, lon,lat, ax, lvls, cmap='Greys', cblbl=r'$E_{??}$', titlestr='', tickintvl=1, norm=None, tick_labels=None, aspect=9, fraction=0.6):

    kwargs_cf = dict(extend='both', norm=norm, cmap=cmap, levels=lvls)
    kwargs_cb = dict(fraction=fraction, aspect=aspect,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])    
    hcf, hcb = sfplt.plotS2field(ax, F, mlon, mlat, kwargs_cf=kwargs_cf, showcb=True, kwargs_cb=kwargs_cb)
    if tick_labels is not None: hcb.set_ticklabels(tick_labels[::tickintvl])
    hcb.set_label(cblbl)
    hcb.ax.xaxis.set_ticks(lvls, minor=True)
#    ax.set_title(titlestr, fontsize=FS, pad=10)
    return hcf

### Generate plots

for deformationType in deformationTypes:


    if DO_DREX:

        blm, F = load_drex_run(deformationType, 1)
        nlm, F = load_drex_run(deformationType, 2)
        vlm, F = load_drex_run(deformationType, 3)
        Nt = F.shape[0]
        Fxz, Fzz = F[:,0,2], F[:,2,2]
        
        blm = blm[:,:nlm_len] # truncated D-Rex states
        nlm = nlm[:,:nlm_len]
        vlm = vlm[:,:nlm_len]
        
        for nn in stepsToPlot:
            print('Selecting D-Rex state for %s at step %i of %i, corresponding to F_zz=%.2f, F_xz=%.2f ...'%(deformationType, nn, Nt, Fzz[nn], Fxz[nn]))
            vlm_est = sf.a4_to_nlm(sf.a4_orth(blm[nn,:], nlm[nn,:]))
            plot_vlm_eval(blm[nn,:],nlm[nn,:],vlm[nn,:],vlm_est, nn, 'DREX-%s-vlm'%(deformationType))
        
        
    if DO_DDM:
    
        dt = 1/100 # Same as used for drex integration in matlab
        
        if deformationType == 'simpleShear':             
            expr = expr_DREX_ss
            T = 0.1225*np.log(2) # time-scale to give same strain per time-step as D-rex was integrated with
            plane = 1
            L = sf.simpleshear_ugrad(plane,T)
            Nt = 3*250
            Fij_t = np.array([sf.simpleshear_F(plane,T, tt*dt) for tt in range(Nt)])
                    
        if deformationType == 'axisymmetricCompression': 
            expr = expr_DREX_uc
            T = 0.115 # time-scale to give same strain per time-step as D-rex was integrated with
            ax, r = 2, 0
            L = sf.pureshear_ugrad(ax,r,T)
            Nt = 2*250
            Fij_t = np.array([sf.pureshear_F(ax,r,T, tt*dt) for tt in range(Nt)])

        Fxz_t, Fzz_t = Fij_t[:,0,2], Fij_t[:,2,2]

        nlm = np.zeros((Nt,nlm_len), dtype=np.complex64)
        blm = nlm.copy()
        vlm = nlm.copy()
        
        fname0 = expr['flist'][0] # initial state 
        fname1 = '%s-%i.csv'%(fname0[:-4],1)
        fname2 = '%s-%i.csv'%(fname0[:-4],2)

        def load_axes(fullfname):
            df = pd.read_csv(fullfname, skiprows=expr['skiprows'], header=None, sep=expr['sep'])
            vj = df.to_numpy()
            return vj

        axes1_0 = load_axes('data/%s/%s'%(expr['path'],fname1))
        axes2_0 = load_axes('data/%s/%s'%(expr['path'],fname2))
        
        mpi = np.zeros((Nt, 3, 3, len(axes1_0))) # m'_i axes (time, m'_i, xyz comp., grain no.) 
        mpi[0,0,:,:] = axes1_0.T
        mpi[0,1,:,:] = axes2_0.T
        mpi[0,2,:,:] = np.cross(axes1_0,axes2_0).T
        
        def mpi_to_nlm(mi):
            r = mi.T
            N,_ = r.shape
            a4 = np.array([ np.einsum('i,j,k,l',r[ii,:],r[ii,:],r[ii,:],r[ii,:], dtype=np.float64) for ii in np.arange(N)]).mean(axis=0) # construct <c^4>
            return sf.a4_to_nlm(a4) 

        ddm_b = sfdsc.DDM(iota=-1, v0=axes1_0)
        ddm_n = sfdsc.DDM(iota=+1, v0=axes2_0)

        Nt = np.amax(stepsToPlot)
        bar = Bar('Integrating DDM for %s until step %i, corresponding to F_zz=%.2f, F_xz=%.2f '%(deformationType, Nt, Fzz_t[Nt], Fxz_t[Nt]), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds')
        for nn in np.arange(1,Nt+1):
            ddm_b.evolve(L, dt)
            ddm_n.evolve(L, dt)
            mpi[nn,0,:,:] = ddm_b.v.T
            mpi[nn,1,:,:] = ddm_n.v.T 
            mpi[nn,2,:,:] = np.cross(ddm_b.v, ddm_n.v).T 
            blm[nn,:] = mpi_to_nlm(mpi[nn,0,:,:])
            nlm[nn,:] = mpi_to_nlm(mpi[nn,1,:,:])
            vlm[nn,:] = mpi_to_nlm(mpi[nn,2,:,:])
            bar.next()
        bar.finish()

        ### Plot selected steps

        shape0 = (len(vlat), len(vlon))

        for nn in stepsToPlot:
    
            print('------------ Plotting step %i ------------'%(nn))

            # DEBUG FOR SINGLE-CRYSTAL BEHAVIOUR (n,b,v are delta functions) ?
            if nn==0: # assume we want to debug
                print('*** DEBUG: OVERRIDING WITH SINGLE CRYSTAL BEHAVIOUR ***')
                x1,x2,x3 = np.eye(3)
                x2,x3 = x3,x2
                if 0:
                    theta_, phi_ = 1*np.pi/4, +0*np.pi/8
                    x1 = sf.rotate_vector(x1, theta_, phi_)
                    x2 = sf.rotate_vector(x2, theta_, phi_)
                    x3 = sf.rotate_vector(x3, theta_, phi_)
                print('b=',x1, 'n=',x2, 'v=',x3)
                mpi = np.zeros((1, 3, 3, 1)) # m'_i axes (time, m'_i, xyz comp., grain no.) 
                mpi[0,0,:,0] = x1
                mpi[0,1,:,0] = x2
                mpi[0,2,:,0] = x3
                blm[nn,:] = mpi_to_nlm(mpi[0,0,:,:])
                nlm[nn,:] = mpi_to_nlm(mpi[0,1,:,:])
                vlm[nn,:] = mpi_to_nlm(mpi[0,2,:,:])

            blm_, nlm_, vlm_ = blm[nn,:], nlm[nn,:], vlm[nn,:] # v = 0 => estimate it from b,n
            mpi_ = mpi[0,:,:,:]


            if PLOT_vlm:

                print('*** Plotting v(r) distr. for DDM vs SDM :: %s, step %i of %i'%(deformationType, nn, Nt))
                ev_v4 = sf.a4_orth(blm_, nlm_) # <v^4> approximation
                vlm_est = sf.a4_to_nlm(ev_v4)
                plot_vlm_eval(blm_,nlm_,vlm_,vlm_est, nn, 'DDM-%s-vlm'%(deformationType))
           

            vlm_ *= 0 # estimimate from blm, nlm below
                
            if PLOT_Vi:

                print('*** Plotting SDM V_i')
                vi_SDM_1d = sf.Vi_elastic_orthotropic(blm_, nlm_, vlm_, alpha,lame_grain,rho, mcolat.flatten(),mlon.flatten())
                vS1_SDM, vS2_SDM, vP_SDM = vi_SDM_1d[0,:].reshape(shape0), vi_SDM_1d[1,:].reshape(shape0), vi_SDM_1d[2,:].reshape(shape0)
                plot_vi_eval(blm_, vP_SDM, vS1_SDM, vS2_SDM, mlon, mlat, nn, 'SDM-%s-vi'%(deformationType))

                print('*** Plotting DDM V_i')
                vi_DDM_1d = sf.Vi_elastic_orthotropic__discrete(mpi_, alpha,lame_grain,rho, mcolat.flatten(),mlon.flatten())
                vS1_DDM, vS2_DDM, vP_DDM = vi_DDM_1d[0,:].reshape(shape0), vi_DDM_1d[1,:].reshape(shape0), vi_DDM_1d[2,:].reshape(shape0)
                plot_vi_eval(blm_, vP_DDM, vS1_DDM, vS2_DDM, mlon, mlat, nn, 'DDM-%s-vi'%(deformationType))
            
            
            if PLOT_Eij:
            
                # r = radial vector, t = theta vector, p = phi vector
                dims = (2, len(vlat), len(vlon))
                Err, Ert, Erp = np.ones(dims), np.ones(dims), np.ones(dims) # r-r longitidinal, r-t shear, r-p shear

                # Discrete (control)
                print('*** Plotting DDM E_ij')
                """
                arr = numpy.zeros((3, 4, 5, 6, 7, 8))
                new_arr = arr.reshape(*arr.shape[:2], -1, *arr.shape[-2:])
                new_arr.shape
                # (3, 4, 30, 7, 8)
                """
                r, t, p = vr.reshape(3,-1), vt.reshape(3,-1), vp.reshape(3,-1)
                tau_rr_, tau_rt_, tau_rp_ = tau_rr.reshape(3,3,-1), tau_rt.reshape(3,3,-1), tau_rp.reshape(3,3,-1)
                Err[0,:,:] = sf.Evw_orthotropic_discrete(mpi_, r, r, tau_rr_, *viscparams_grain).reshape(shape0)
                Ert[0,:,:] = sf.Evw_orthotropic_discrete(mpi_, r, t, tau_rt_, *viscparams_grain).reshape(shape0)
                Erp[0,:,:] = sf.Evw_orthotropic_discrete(mpi_, r, p, tau_rp_, *viscparams_grain).reshape(shape0)

                print('*** Plotting SDM E_ij')
                for ii,_ in enumerate(vlat):
                    for jj,_ in enumerate(vlon):
                    
                        r,t,p = vr[:,ii,jj], vt[:,ii,jj], vp[:,ii,jj]
                            
                        # Spectral estimate
                        Err[1,ii,jj] = sf.Evw_orthotropic(blm_, nlm_, vlm_,  r, r, tau_rr[:,:,ii,jj], *viscparams_grain)
                        Ert[1,ii,jj] = sf.Evw_orthotropic(blm_, nlm_, vlm_,  r, t, tau_rt[:,:,ii,jj], *viscparams_grain) 
                        Erp[1,ii,jj] = sf.Evw_orthotropic(blm_, nlm_, vlm_,  r, p, tau_rp[:,:,ii,jj], *viscparams_grain)
                
                plot_Eij_eval(blm_, nlm_, vlm_, Err[0,:,:],Ert[0,:,:],Erp[0,:,:], mlon, mlat, nn, 'DDM-%s-Eij'%(deformationType))
                plot_Eij_eval(blm_, nlm_, vlm_, Err[1,:,:],Ert[1,:,:],Erp[1,:,:], mlon, mlat, nn, 'SDM-%s-Eij'%(deformationType))

#                    a2v_nlm, a2v_mi, a4v_nlm, a4v_mi = sf.orthstructs(blm_, nlm_, vlm_, mi)
#                    
#                    # *** IS THIS THE SOLUTION?? k=2 (v moments) show large discrepancy! ***
#                    kk = 2 
#                    print('a2v(%i,:) diff'%(kk))
#                    print(a2v_mi[kk,:]-a2v_nlm[kk,:])
#                    print('a4v(%i,:,:) diff'%(kk))
#                    print(a4v_mi[kk,:,:]-a4v_nlm[kk,:,:])
#                    
#                    code.interact(local=locals())                    
                           
