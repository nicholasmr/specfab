#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

import copy, sys, os, code # code.interact(local=locals())

import numpy as np
from scipy import interpolate
import scipy.special as sp

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import common as sfcom
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, rc, patches

#colorbg = "#f7f6ee"
os.system('mkdir -p frames')    
exp = str(sys.argv[1])

### Determine fabric evolution

L = 18 # Spectral truncation
lm, nlm_len = sf.init(L)

# Model LROT only
kwargs_CPOevo = dict(iota=+1, Gamma0=None, Lambda=None, nu=1)

# uniaxial compression
strain_target = -0.92
Nt_uc = 100
nlm_uc, F_uc, time, ugrad = sfint.lagrangianparcel(sf, dict(type='ps', axis=0, T=+1, r=0), strain_target, Nt=Nt_uc, **kwargs_CPOevo) 
strain_uc = np.array([sf.F_to_strain(F_uc[tt,:,:])[0,0] for tt in range(Nt_uc+1)])

# uniaxial extension
strain_target = +2
Nt_ue = 100
nlm_ue, F_ue, time, ugrad = sfint.lagrangianparcel(sf, dict(type='ps', axis=0, T=-1, r=0), strain_target, Nt=Nt_ue, **kwargs_CPOevo) 
strain_ue = np.array([sf.F_to_strain(F_ue[tt,:,:])[0,0] for tt in range(Nt_uc+1)])

### Determine enhancement factors and eigenvalues

grain_params = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)

Eij_ue = sf.Eij_tranisotropic_arr(nlm_ue, *sfcom.xi_tile(Nt_ue+1), *grain_params) 
_, eigvals_ue = sfcom.eigenframe(nlm_ue)

Eij_uc = sf.Eij_tranisotropic_arr(nlm_uc, *sfcom.xi_tile(Nt_uc+1), *grain_params) 
_, eigvals_uc = sfcom.eigenframe(nlm_uc)

### Plot

Nt      = Nt_uc      if exp=='uc' else Nt_ue 
nlm     = nlm_uc     if exp=='uc' else nlm_ue 
F       = F_uc       if exp=='uc' else F_ue
eigvals = eigvals_uc if exp=='uc' else eigvals_ue
strain  = strain_uc  if exp=='uc' else strain_ue

for tt in range(Nt+1):
    
    scale = 0.4
    fig = plt.figure(figsize=(13*scale,10*scale))
    
    gs_master = gridspec.GridSpec(1, 2, width_ratios=[0.8,1])
    gs_master.update(top=0.91, bottom=0.13, left=0.095, right=1-0.05, hspace=-0.05, wspace=-0.075)
    
    gsleft = gs_master[0,0].subgridspec(2, 1,  height_ratios=[1,0.3]) 
    axE = fig.add_subplot(gsleft[0, 0])
    
    gs = gs_master[0,1].subgridspec(2, 1,  height_ratios=[1,0.8], hspace=-0.05)
    axParcel = fig.add_subplot(gs[0, 0], projection='3d')
    
    geo, prj = sfplt.getprojection(rotation=50+180, inclination=50)
          
    axODF = fig.add_subplot(gs[1, 0], projection=prj)
    axODF.set_global()
    
    ### Plot enhancement factors

    xlims = [-1,2]
    ylims = [0,2.5]            

    axE.add_patch(plt.Rectangle((xlims[0],1), np.diff(xlims)[0], ylims[1]-1, color='#deebf7'))
    axE.add_patch(plt.Rectangle((xlims[0],0), np.diff(xlims)[0], 1, color='#fee0d2'))

    lw = 1.5
    axE.plot(strain[[tt,tt]], ylims, '-', c='0.5', lw=lw)

    # Eij = (Exx,Eyy,Ezz,Eyz,Exz,Exy)
    axE.plot(strain_ue, Eij_ue[:,0], '--k', lw=lw, label=r'$E_{xx}$')
    axE.plot(strain_ue, Eij_ue[:,2],  ':k', lw=lw, label=r'$E_{zz}$')
    axE.plot(strain_ue, Eij_ue[:,4], '-k',  lw=lw, label=r'$E_{xz}$')
    axE.plot(strain_ue, Eij_ue[:,3], '-.k', lw=lw, label=r'$E_{yz}$')

    axE.plot(strain_uc, Eij_uc[:,0], '--k', lw=lw)
    axE.plot(strain_uc, Eij_uc[:,2],  ':k', lw=lw)
    axE.plot(strain_uc, Eij_uc[:,4], '-k',  lw=lw)
    axE.plot(strain_uc, Eij_uc[:,3], '-.k', lw=lw)
    
    axE.text(0.02, 1.50, '{\\bf Softer than}\n{\\bf isotropy}', fontsize=FS-2.5, color='#08519c', ha='center', ma='center')
    axE.text(0.02, 0.20, '{\\bf Harder than}\n{\\bf isotropy}', fontsize=FS-2.5, color='#a50f15', ha='center', ma='center')

    dx=1
    axE.set_xticks(np.arange(xlims[0],xlims[1]+dx, dx))
    axE.set_xticks(np.arange(xlims[0],xlims[1]+dx, dx/2), minor=True)
    axE.set_xlim(xlims)

    dy=1
    axE.set_yticks(np.arange(ylims[0],ylims[1]+dy,dy))
    axE.set_yticks(np.arange(ylims[0],ylims[1]+dy,dy/4), minor=True)
    axE.set_ylim(ylims)
    
    axE.set_xlabel(r'$\epsilon_{xx}$')
    axE.set_ylabel('$E_{ij}$')
    axE.set_title('Enhancement factors', pad=8, fontsize=FS-0.5)
    #
    legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'handlelength':1.9, 'labelspacing':0.2}
    hleg = axE.legend(loc=1, fontsize=FS-0.5, **legkwargs)
    hleg.get_frame().set_linewidth(0.7);
    
    ### Plot annotations
    
    x0, y0, dy = -0.125, -0.38, 0.09
    axE.text(x0, y0, r'{\bf specfab}', transform=axE.transAxes, fontsize=FS, ha='left')
    axE.text(x0, y0-1*dy, r'github.com/nicholasmr/specfab', transform=axE.transAxes, fontsize=FS-1.5, ha='left')
    axE.text(x0, y0-2.0*dy, r'Lattice rotation demo', transform=axE.transAxes, fontsize=FS-1.5, ha='left')
    
    x0, y0, dy = +1.875, -0.15, 0.08
    axE.text(x0,y0+2*dy, r'$\lambda_1 = %.2f$'%(eigvals[tt,0]), transform=axE.transAxes, fontsize=FS-1, ha='left')
    axE.text(x0,y0+1*dy, r'$\lambda_2 = %.2f$'%(eigvals[tt,1]), transform=axE.transAxes, fontsize=FS-1, ha='left')
    axE.text(x0,y0+0*dy, r'$\lambda_3 = %.2f$'%(eigvals[tt,2]), transform=axE.transAxes, fontsize=FS-1, ha='left')
    
    ### Plot parcel geometry
    
    sfplt.plotparcel(axParcel, F[tt,:,:], azim=25, scale=1, axscale=1.8) 
    
    ### Plot ODF

    sfplt.plotODF(nlm[tt,:], lm, axODF, lvlset=[np.linspace(0,0.4,9), lambda x,p:'%.1f'%x], cblabel='ODF')
    sfplt.plotcoordaxes(axODF, geo, axislabels='vuxi', color='k')
    
    ### Save frame in correct order (file name number) for generating animation below

    fout = 'frames/frame-%s-%i.png'%(exp,tt)
    print('Saving %i :: %s'%(tt,fout))    
    plt.savefig(fout, transparent=True, dpi=300)
    plt.close()

