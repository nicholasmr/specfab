# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2024

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
FS = sfplt.setfont_tex(fontsize=11)
FSLEG = FS-1

L = 4
lm, nlm_len = sf.init(L)

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

DO_DREX = 1 # D-Rex runs (to be re-loaded)
DO_DDM  = 1 # Discrete Directors Method 

PLOT_vlm = 0
PLOT_Vi  = 1

# Experiments to process
deformationTypes = ['simpleShear','axisymmetricCompression']
#deformationTypes = ['simpleShear',]
#deformationTypes = ['axisymmetricCompression',]

stepsToPlot = [15,] # only 15 is used for papar plots
#stepsToPlot = [0,] # debug, single crystal behaviour

# r = radial vector, t = theta vector, p = phi vector
(vr,vt,vp, rr,rt,rp, tau_rr,tau_rt,tau_rp, vlon,vlat, mlon,mlat) = sfdsc.idealstressstate(lonres=80, latres=40)
mcolat = sfdsc.lat2colat(mlat)

### Functions

lvlset = [np.linspace(0,0.2,9), lambda x,p:'%.1f'%x] # [np.linspace(0.0, 0.2, 5), lambda x,p:'%.1f'%x]

def plot_vlm_eval(blm,nlm,vlm,vlm_est, ii, deformationType):
   
    size = 2.0
    fig = plt.figure(figsize=(1.8*size,0.7*size))
    gs = fig.add_gridspec(ncols=5, nrows=1, width_ratios=[1,1,1,0.15,1], \
                    wspace=0.27, top=1, bottom=0.325, left=0.02, right=1-0.02)
    axb, axn, ax1, ax2 = [fig.add_subplot(gs[0,ii], projection=prj) for ii in [0,1,2,4]] 

    kwargs = dict(lvlset=lvlset, cbtickintvl=4, cbpad=0.1, cbaspect=8, nchunk=0)
    sfplt.plotODF(blm, lm, axb, cblabel=r'$b/N$', cmap='Reds', **kwargs)
    sfplt.plotODF(nlm, lm, axn, cblabel=r'$n/N$', cmap='Blues', **kwargs)
    sfplt.plotODF(vlm, lm, ax1, cblabel=r'$v/N$', cmap='Greens', **kwargs)
    sfplt.plotODF(vlm_est, lm, ax2, cblabel=r'$v/N$', cmap='Greens', **kwargs)

    for ax in (axb,axn,ax1,ax2):
        ax.set_global()
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', color='k', fontsize=FS)
    
    fname = 'texplots/%s-%i.pdf'%(deformationType, ii)
    plt.savefig(fname, transparent=True)
    os.system('pdfcrop %s %s'%(fname,fname))

def plot_vi_eval(blm, vP,vS1,vS2, lon, lat, ii, deformationType):

    size = 2.0
    fig = plt.figure(figsize=(3/4*2.3*size,0.95*size))
    gs = fig.add_gridspec(ncols=3, nrows=1)
    ax1, ax2, ax3 = [fig.add_subplot(gs[0,ii], projection=prj) for ii in range(3)] 

    lvlsP = np.linspace(7.4, 9.2, 7)
    lvlsS1 = np.linspace(4.6, 5.2, 7) # np.linspace(4.7, 5.3, 9)
    lvlsS2 = lvlsS1
    plot_vi(vP*1e-3,  lon,lat, ax=ax1, cblabel='$V_{\mathrm{P}}$ (\SI{}{\kilo\metre\per\second})',  lvls=lvlsP)
    plot_vi(vS1*1e-3, lon,lat, ax=ax2, cblabel='$V_{\mathrm{S1}}$ (\SI{}{\kilo\metre\per\second})', lvls=lvlsS1)
    plot_vi(vS2*1e-3, lon,lat, ax=ax3, cblabel='$V_{\mathrm{S2}}$ (\SI{}{\kilo\metre\per\second})', lvls=lvlsS2)

    for ax in (ax1,ax2,ax3):
        ax.set_global()
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', color='k', fontsize=FS)
    
    plt.tight_layout()
    fname = 'texplots/%s-%i.pdf'%(deformationType, ii)
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
            ugrad = sf.simpleshear_ugrad(plane,T)
            Nt = 3*250
            Fij_t = np.array([sf.simpleshear_F(plane,T, tt*dt) for tt in range(Nt)])
                    
        if deformationType == 'axisymmetricCompression': 
            expr = expr_DREX_uc
            T = 0.115 # time-scale to give same strain per time-step as D-rex was integrated with
            ax, r = 2, 0
            ugrad = sf.pureshear_ugrad(ax,r,T)
            Nt = 2*250
            Fij_t = np.array([sf.pureshear_F(ax,r,T, tt*dt) for tt in range(Nt)])

        Fxz_t, Fzz_t = Fij_t[:,0,2], Fij_t[:,2,2]
        
        fname0 = expr['flist'][0] # initial state 
        fname1 = '%s-%i.csv'%(fname0[:-4],1)
        fname2 = '%s-%i.csv'%(fname0[:-4],2)

        def load_axes(fullfname):
            df = pd.read_csv(fullfname, skiprows=expr['skiprows'], header=None, sep=expr['sep'])
            return df.to_numpy() # vj
            
        axes1_0 = load_axes('data/%s/%s'%(expr['path'],fname1))
        axes2_0 = load_axes('data/%s/%s'%(expr['path'],fname2))

        # All-in-one Euler integration
        b0 = axes1_0 # grain no., xyz
        n0 = axes2_0
        N = np.amax(stepsToPlot) + 1
        eps = np.array([(ugrad+ugrad.T)/2]*N)
        omg = np.array([(ugrad-ugrad.T)/2]*N)
        bi = sf.ri_LROT(b0, dt, N, eps,omg, -1) # time, grain no., xyz
        ni = sf.ri_LROT(n0, dt, N, eps,omg, +1)
        vi = np.array([np.cross(bi[nn,:,:], ni[nn,:,:]) for nn in range(N)])

        Ngrains = bi.shape[1] # no grains
        wi = np.ones(Ngrains)/Ngrains # equal weights
        blm = np.array([sf.ri_to_nlm(bi[nn], wi, L) for nn in range(N)])
        nlm = np.array([sf.ri_to_nlm(ni[nn], wi, L) for nn in range(N)])
        vlm = np.array([sf.ri_to_nlm(vi[nn], wi, L) for nn in range(N)])

        ### Plot selected steps

        shape0 = (len(vlat), len(vlon))

        for nn in stepsToPlot:
    
            print('------------ Plotting step %i ------------'%(nn))

            # DEBUG FOR SINGLE-CRYSTAL BEHAVIOUR (n,b,v are delta functions) ?
            if nn==0: # assume we want to debug
                print('*** DEBUG: OVERRIDING WITH SINGLE CRYSTAL BEHAVIOUR ***')
                x1,x3,x2 = np.eye(3)
                if 0:
                    theta_, phi_ = 1*np.pi/4, +0*np.pi/8
                    x1 = sf.rotate_vector(x1, theta_, phi_)
                    x2 = sf.rotate_vector(x2, theta_, phi_)
                    x3 = sf.rotate_vector(x3, theta_, phi_)
                print('b=',x1, 'n=',x2, 'v=',x3)
                bi, ni, vi = x1, x2, x3
                blm = np.array([sf.ri_to_nlm(np.array(bi,ndmin=2), wi, L),])
                nlm = np.array([sf.ri_to_nlm(np.array(ni,ndmin=2), wi, L),])
                vlm = np.array([sf.ri_to_nlm(np.array(vi,ndmin=2), wi, L),])
                
            blm_, nlm_, vlm_ = blm[nn,:], nlm[nn,:], vlm[nn,:] # v=0 => estimate it from b,n

            if PLOT_vlm:

                print('*** Plotting v(r) distr. for DDM vs SDM :: %s, step %i of %i'%(deformationType, nn, Nt))
                ev_v4 = sf.a4_orth(blm_, nlm_) # <v^4> approximation
                vlm_est = sf.a4_to_nlm(ev_v4)
                plot_vlm_eval(blm_,nlm_,vlm_,vlm_est, nn, 'DDM-%s-vlm'%(deformationType))
           

            if PLOT_Vi:

                vlm_ *= 0 # estimimate from (blm,nlm)

                print('*** Plotting SDM V_i')
                vi_SDM_1d = sf.Vi_elastic_orthotropic(blm_, nlm_, vlm_, alpha,lame_grain,rho, mcolat.flatten(),mlon.flatten())
                vS1_SDM, vS2_SDM, vP_SDM = vi_SDM_1d[0,:].reshape(shape0), vi_SDM_1d[1,:].reshape(shape0), vi_SDM_1d[2,:].reshape(shape0)
                plot_vi_eval(blm_, vP_SDM, vS1_SDM, vS2_SDM, mlon, mlat, nn, 'SDM-%s-vi'%(deformationType))

                print('*** Plotting DDM V_i')
                mpi = np.array([bi[nn].T,ni[nn].T,vi[nn].T])
                vi_DDM_1d = sf.Vi_elastic_orthotropic__discrete(mpi, alpha,lame_grain,rho, mcolat.flatten(),mlon.flatten())
                vS1_DDM, vS2_DDM, vP_DDM = vi_DDM_1d[0,:].reshape(shape0), vi_DDM_1d[1,:].reshape(shape0), vi_DDM_1d[2,:].reshape(shape0)
                plot_vi_eval(blm_, vP_DDM, vS1_DDM, vS2_DDM, mlon, mlat, nn, 'DDM-%s-vi'%(deformationType))              
                           
