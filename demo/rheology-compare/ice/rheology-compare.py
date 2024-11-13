# N. M. Rathmann <rathmann@nbi.ku.dk>

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import common as sfcom
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSSMALL = FS-1

#----------------------
# Parameters
#----------------------

# Bulk flow law
nglen = 3
Aglen = 3.5e-26 # A(T=-25 deg.)

#----------------------
# Experiment
#----------------------

T_EXP_CC = 1 # confined vertical compression 
T_EXP_SS = 2 # vertical shear

# Select experiment
T_EXP = T_EXP_CC
T_EXP = T_EXP_SS

#----------------------
# Experiment definitions
#----------------------

Nt = 200
L  = 18

if T_EXP==T_EXP_SS:
    i,j = 0,2 # components of interest
    mod = dict(type='ss', plane=1)
    strain_target = np.deg2rad(80.5) # simulate  parcel deformation until this target strain
    strain_thres = 2.5 # assume DDRX dominates when reached this strain threshold
    ODF_strains = [0, 1, 2, 3] # show ODFs at these strains
    xlims = [0, 3.0]
    ylims = [0.7, 2.7]
                
if T_EXP==T_EXP_CC:
    i = j = 2 # components of interest
    mod = dict(type='ps', r=+1, axis=i)
    strain_target = -0.97 # simulate  parcel deformation until this target strain
    strain_thres = -0.7 # assume DDRX dominates when reached this strain threshold
    ODF_strains = [0, -0.3, -0.5, -0.9] # show ODFs at these strains
    xlims = [0, -1]
    ylims = [0.3, 3]
        
#----------------------
# CPO evolution
#----------------------

lm, nlm_len = sf.init(L)
nlm = np.zeros((Nt+1,nlm_len), dtype=np.complex64) # The expansion coefficients

kwargs_LROT = dict(iota=1, Gamma0=None, nu=1) #      
nlm[:,:], F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, **kwargs_LROT)

# Strain history
strain = np.array([sf.F_to_strain(F[tt,:,:])[i,j] for tt in range(Nt+1)])
D, W = sf.ugrad_to_D_and_W(ugrad)
S = D/np.linalg.norm(D) # stress coaxial with strain-rate

# Model DDRX from strain = strain_thres until strain_target
# (this is done by splicing DDRX simulation results into array containing LROT simulation results)
I = np.argmin(np.abs(strain-strain_thres)) # starting state
kwargs_DDRX = dict(iota=None, Gamma0=10, nu=1)
nlm_, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, nlm0=nlm[I,:], **kwargs_DDRX)
nlm[I:,:] = nlm_[:(Nt-I+1),:]

#----------------------
# Determine eigenenhancements etc.
#----------------------

# G=Glen, R=Rathmann & Lilien, P=Pettit, M=Martin
Y_G = sf.rheo_fwd_isotropic(S, Aglen,nglen)[i,j]
dims = (Nt+1)
Y_R, Y_P, Y_M = np.zeros(dims), np.zeros(dims), np.zeros(dims)

grain_params = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
ei, eig = sfcom.eigenframe(nlm)
e1,e2,e3 = ei[:,:,0], ei[:,:,1], ei[:,:,2]
Eij = sf.Eij_tranisotropic_arr(nlm, e1,e2,e3, *grain_params) # eigenenhancements
    
# Euler integration 
for tt in np.arange(0,Nt+1):
    # Y for fabric at constant strain rate
    args = (S,Aglen,nglen, e1[tt],e2[tt],e3[tt], Eij[tt])
    Y_R[tt] = sf.rheo_fwd_orthotropic(*args)[i,j]
    Y_P[tt] = sf.rheo_fwd_orthotropic_Pettit(*args)[i,j]
    Y_M[tt] = sf.rheo_fwd_orthotropic_Martin(*args)[i,j]

#----------------------
# Plot
#----------------------

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec

geo, prj = sfplt.getprojection(rotation=50+180, inclination=50)
ODF_tsteps = np.array([np.argmin(np.abs(strain-thres)) for thres in ODF_strains])

### Setup figure

scale = 0.87
fig = plt.figure(figsize=(7.0*scale,5*scale), constrained_layout=False)
gs = gridspec.GridSpec(2, 2, height_ratios=[1,0.55])
gs.update(hspace=0.53, wspace=0.26, left=0.10, right=0.97, top=0.98, bottom=0.12)
ax_Y = fig.add_subplot(gs[0, :])
gs_parcels = gridspec.GridSpecFromSubplotSpec(1, len(ODF_tsteps), subplot_spec=gs[1,:], width_ratios=[1,1,1,1], hspace=0.6, wspace=+0.2)
ax_ODF = [fig.add_subplot(gs_parcels[0, ii], projection=prj) for ii,nn in enumerate(ODF_tsteps)]

### Plot relative strain-rates

# Background patches

ax_Y.add_patch(plt.Rectangle((strain[0],0), strain_thres,  10, color=sfplt.c_vlblue))
ax_Y.add_patch(plt.Rectangle((strain_thres,0), strain[-1], 10, color=sfplt.c_vlred))

# Glen's line

cglen = '0.35'
if T_EXP==T_EXP_SS: xglen,yglen = 1.95, 1.05
if T_EXP==T_EXP_CC: xglen,yglen = -0.625, 1.08
ax_Y.plot(strain, strain*0+1, '-', lw=1.2,  color=cglen)
ax_Y.text(xglen, yglen, "Glen's law", color=cglen, ha='center', fontsize=FSSMALL)

# Model relative lines

Yn_R, Yn_P, Yn_M = Y_R/Y_G, Y_P/Y_G, Y_M/Y_G
ax_Y.plot(strain, Yn_R, '-',  c='k', label=r"Unapprox.")
ax_Y.plot(strain, Yn_M, '--', c='k', label=r"M09")
ax_Y.plot(strain, Yn_P, ':',  c='k', label=r"P07")

if 0:
    delta_Martin = 100*np.divide(np.abs(Yn_M-Yn_R),Yn_R)
    delta_Pettit = 100*np.divide(np.abs(Yn_P-Yn_R),Yn_R)
    print('max(Martin-True)=%.2f'%(np.amax(delta_Martin)))
    print('max(Pettit-True)=%.2f'%(np.amax(delta_Pettit)))
    print('(Martin-True)_{-1}=%.2f'%(delta_Martin[-1]))
    print('(Pettit-True)_{-1}=%.2f'%(delta_Pettit[-1]))

# Axis labels

ax_Y.set_xlabel(r'$\epsilon$')
ijstr = '%s%s'%('x' if i==0 else 'z', 'x' if j==0 else 'z')
ax_Y.set_ylabel(r'$\dot{\epsilon}_{%s}/\dot{\epsilon}^{\mathrm{Glen}}_{%s}$'%(ijstr,ijstr))

# Legend

legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':1, 'ncol':1, 'handlelength':1.34, 'labelspacing':0.25, 'fontsize':FSSMALL}
leg = ax_Y.legend(loc=4, **legkwargs) # bbox_to_anchor=(1,0.88)
leg.get_frame().set_linewidth(0.8);
sfplt.panellabel(ax_Y, 2, r'{\bf a}')

# Axis ticks and limits

if T_EXP==T_EXP_SS: xticks = np.arange(0,20,0.25) 
if T_EXP==T_EXP_CC: xticks = np.arange(-1,0,0.1)
ax_Y.set_xticks(xticks, minor=True)  
ax_Y.set_xlim(xlims)

ax_Y.set_yticks(np.arange(0,10+1e-3,0.5))
ax_Y.set_yticks(np.arange(0,10+1e-3,0.1), minor=True)
ax_Y.set_ylim(ylims)

# Set bg patches identifies

dx = 0.016
y0 = ylims[-1]-0.075*abs(np.diff(ylims)[0])
ax_Y.text(strain_thres*(1-dx), y0, r'{\bf Lattice rotation}', c=sfplt.c_dblue, ha='right', va='center', fontsize=FSSMALL)
ax_Y.text(strain_thres*(1+dx), y0, r'{\bf DDRX}',             c=sfplt.c_dred,  ha='left',  va='center', fontsize=FSSMALL)

### Plot ODFs

for ii,nn in enumerate(ODF_tsteps):

    ax = ax_ODF[ii]
    ax.set_global()
    sfplt.plotODF(nlm[nn,:], lm, ax, cblabel=r'ODF', lvlset=[np.linspace(0,0.41,9), lambda x,p:'%.1f'%x])
    ax.set_title(r'$\epsilon=%.1f$'%(strain[nn]), fontsize=FS)
    sfplt.panellabel(ax,2,r'{\bf %s}'%(chr(ord('b')+ii)), bbox=(-0.25,1.35))

    if ii==0: 
        sfplt.plotcoordaxes(ax, geo, axislabels='vuxi', color='k')        
    elif ii>0:    
#        w,v = np.linalg.eig(tau_RL[nn,:,:]) # debug: principal stress directions
#        for ei in (v[:,0],v[:,1],v[:,2]): # debug: principal stress directions
        for kk, ei in enumerate((e1[nn,:],e2[nn,:],e3[nn,:])):
            ck = [sfplt.c_dred, sfplt.c_dgreen, sfplt.c_dblue]
            kwargs = dict(marker='.', ms=7, markerfacecolor=ck[kk], markeredgecolor=ck[kk], markeredgewidth=1.0, transform=geo)
            sfplt.plotS2point(ax, +ei, **kwargs)
            sfplt.plotS2point(ax, -ei, **kwargs)

### Plot parcel deformation as inset

if T_EXP==T_EXP_SS: parcel_tsteps = ODF_tsteps[[0,1]]
if T_EXP==T_EXP_CC: parcel_tsteps = ODF_tsteps[[0,2]]

for ii,nn in enumerate(parcel_tsteps):

    # Parcel locations (normalized figure coords, so must be set manually like this)
    if T_EXP==T_EXP_SS: 
        y0 = 0.7
        pc = np.array([[0.106,y0],[0.28,y0]])
        
    if T_EXP==T_EXP_CC: 
        y0 = 0.7
        pc = np.array([[0.106,y0],[0.45,y0]])

    axs = 0.14
    ax_sub = fig.add_axes([pc[ii,0], pc[ii,1], axs, axs], projection='3d')
    ax_sub.patch.set_alpha(0.0)
    lw = 0.7
    sfplt.plotparcel(ax_sub, F[nn,:,:], lw=lw, lwinit=lw, fonttex=True)
    ax_sub.set_title(r'$\epsilon=%.1f$'%(strain[nn]), fontsize=FSSMALL,  x=0.5, y=0.92)
        
### Save figure

T_EXP_STR = {T_EXP_CC:'cc', T_EXP_SS:'ss'}
fname = 'rheology-compare-%s.png'%(T_EXP_STR[T_EXP])
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=175)

