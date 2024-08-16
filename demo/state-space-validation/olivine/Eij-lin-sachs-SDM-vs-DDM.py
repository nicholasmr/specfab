# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp
from progress.bar import Bar

sys.path.insert(0, '..')
from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import discrete as sfdsc
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSSMALL = FS-1

ENABLE_DDM = 1

SHOW_ALT_EIJ_COMPONENTS = False

#----------------------
# Experiment
#----------------------

T_EXP_UC = 1 # confined vertical compression 
T_EXP_SS = 2 # vertical shear

# Select experiment
T_EXP = T_EXP_UC
#T_EXP = T_EXP_SS

#----------------------
# Experiment definitions
#----------------------

Nt = 300
L  = 10

if T_EXP==T_EXP_SS:
    i,j = 0,2 # Fij components of interest
    mod = dict(type='ss', plane=1)
    strain_target = np.deg2rad(79) # simulate  parcel deformation until this target strain
    ODF_strains = [1, 4] # show ODFs at these strains
    xlims = [0, 5]
    ylims = [0, 3.1]
    expr = expr_DREX_ss # for DDM init state
                
if T_EXP==T_EXP_UC:
    i = j = 2 # Fij components of interest
    mod = dict(type='ps', r=0, axis=i)
    strain_target = -0.99 # simulate  parcel deformation until this target strain
    ODF_strains = [0.2, 0.5] # show ODFs at these strains
    xlims = [-1, 0]
    ylims = [0, 2.5]
    expr = expr_DREX_uc # for DDM init state
        
def strainax(F):
    if T_EXP==T_EXP_SS: return f_gammaxz(F)
    if T_EXP==T_EXP_UC: return f_strainzz(F)
    
#----------------------
# CPO evolution
#----------------------

### SDM

lm, nlm_len = sf.init(L)
blm = np.zeros((Nt+1,nlm_len), dtype=np.complex64) 
nlm = np.zeros((Nt+1,nlm_len), dtype=np.complex64) 

# Slip plane normal distribution
kwargs_LROT = dict(iota=+1, Gamma0=None, nu=1)
nlm[:,:], F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, **kwargs_LROT)

# Slip direction distribution
kwargs_LROT = dict(iota=-1, Gamma0=None, nu=1)
blm[:,:], F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, **kwargs_LROT)

### DDM
            
fname0 = expr['flist'][0] # initial state from DREX ensemble
fname1 = '%s-%i.csv'%(fname0[:-4],1)
fname2 = '%s-%i.csv'%(fname0[:-4],2)

def load_axes(fullfname):
    df = pd.read_csv(fullfname, skiprows=expr['skiprows'], header=None, sep=expr['sep'])
    vj = df.to_numpy()
    return vj

axes1_0 = load_axes('data/%s/%s'%(expr['path'],fname1))
axes2_0 = load_axes('data/%s/%s'%(expr['path'],fname2))

mpi = np.zeros((Nt+1, 3, 3, len(axes1_0))) # m'_i axes (time, m'_i, xyz comp., grain no.) 
mpi[0,0,:,:] = axes1_0.T
mpi[0,1,:,:] = axes2_0.T
mpi[0,2,:,:] = np.cross(axes1_0,axes2_0).T

ddm_b = sfdsc.DDM(iota=-1, v0=axes1_0)
ddm_n = sfdsc.DDM(iota=+1, v0=axes2_0)

if ENABLE_DDM:

    dt = time[1]-time[0]
    nn = 0
    bar = Bar('DDM :: Nt=%i :: dt=%.2e '%(Nt,dt), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds')
    while nn < Nt:
        nn += 1
        ddm_b.evolve(ugrad, dt)
        ddm_n.evolve(ugrad, dt)
        
        mpi[nn,0,:,:] = ddm_b.v.T
        mpi[nn,1,:,:] = ddm_n.v.T 
        mpi[nn,2,:,:] = np.cross(ddm_b.v, ddm_n.v).T 

        bar.next()
    bar.finish()
        
#----------------------
# Determine eigenenhancements etc.
#----------------------

dims = (Nt+1,3)
e1,e2,e3,eig = np.zeros(dims),np.zeros(dims),np.zeros(dims),np.zeros(dims)
Eij = np.zeros((Nt+1,6))
Eij_dsc = np.zeros((Nt+1,6))

ijpars = [(0,0), (1,1), (2,2), (1,2), (0,2), (0,1)] # Voigt ordering

grain_params = sfconst.olivine['viscoplastic']['linear'] # Optimal n'=1 Sachs grain parameters
print('Grain params: ', grain_params)

bar = Bar('Eij :: '%(), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds')
for nn in np.arange(0,Nt+1):

    ### SDM 
    
    nlm_ = nlm[nn,:] 
    blm_ = blm[nn,:] 
    vlm_ = 0*nlm_ # infer from nlm, blm
    
    e1[nn,:], e2[nn,:], e3[nn,:], eig[nn,:] = sf.frame(nlm_, 'e') # eigenframe
    mi = np.zeros((3,3))
    mi[0,:] = e1[nn,:]
    mi[1,:] = e2[nn,:]
    mi[2,:] = e3[nn,:]

    nlm_1, nlm_2, nlm_3 = blm_, nlm_, vlm_ # A-type olivine
    
    for kk, (i_,j_) in enumerate(ijpars):
        mi_, mj_ = mi[i_,:],mi[j_,:]
        mimj = np.tensordot(mi_, mj_, axes=0)
        if i_==j_: tauhat = np.identity(3) - 3*mimj
        else:    tauhat = mimj + mimj.T
        Eij[nn,kk] = sf.Evw_orthotropic(nlm_1,nlm_2,nlm_3, mi_, mj_,tauhat, *grain_params)
        if ENABLE_DDM: Eij_dsc[nn,kk] = sf.Evw_orthotropic_discrete(mpi[nn,:,:,:], mi_, mj_,tauhat, *grain_params)

    bar.next()
bar.finish()

#----------------------
# Plot
#----------------------

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec

geo, prj = sfplt.getprojection(rotation=50+180, inclination=50)
ODF_tsteps = np.array([np.argmin(np.abs(F[:,i,j]-thres)) for thres in ODF_strains])

### Setup figure

scale = 0.87
fig = plt.figure(figsize=(7.0*scale,4.35*scale), constrained_layout=False)
gs = gridspec.GridSpec(2, 2, height_ratios=[1,0.475])
gs.update(hspace=0.53, wspace=0.26, left=0.11, right=0.98, top=0.96, bottom=0.12)
ax_Y = fig.add_subplot(gs[0, :])
gs_parcels = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[1,:], width_ratios=[1,1,1.5,1,1], hspace=0.6, wspace=+0.2)
ax_ODF = [fig.add_subplot(gs_parcels[0, ii], projection=prj) for ii in [0,1, 3,4]]

### Plot relative strain-rates

# Isotropic line

ciso = '0.15'
if T_EXP==T_EXP_SS: xiso,dyiso = 2.5, 0.2
if T_EXP==T_EXP_UC: xiso,dyiso = 0.03-1, 0.16
ax_Y.plot(xlims, [1,1], ':', lw=1.2,  color=ciso)
#ax_Y.text(xiso, yiso, "Isotropic olivine", color=ciso, ha='center', va='center', fontsize=FSSMALL)

#X = xlims
#ax_Y.add_patch(plt.Rectangle((X[0],1), X[-1], ylims[-1], color='#d9f0d3'))
#ax_Y.add_patch(plt.Rectangle((X[0],1), X[-1], -1,        color='#e7d4e8'))

ax_Y.text(xiso, 1+dyiso, r'{$\uparrow$   Softer than isotropy}', c=ciso, ha='left', va='center', fontsize=FSSMALL)
ax_Y.text(xiso, 1-dyiso, r'{$\downarrow$ Harder than isotropy}', c=ciso, ha='left', va='center', fontsize=FSSMALL)

# Model relative lines
I11, I22, I33 = 0,1,2
I23, I13, I12 = 3,4,5

ax_Y.plot(strainax(F), Eij[:,I13], '-',  c='k', label=r"$E_{13}$ SDM")
ax_Y.plot(strainax(F), Eij[:,I11], '--', c='k', label=r"$E_{11}$ SDM")

if SHOW_ALT_EIJ_COMPONENTS:
    ax_Y.plot(strainax(F), Eij[:,I22], '-.', c='0.5',  label=r"$E_{22}$ SDM")
    ax_Y.plot(strainax(F), Eij[:,I33], ':',  c='0.75', label=r"$E_{33}$ SDM")

ax_Y.plot(strainax(F), Eij_dsc[:,I13], '-',  c='tab:green', label=r"$E_{13}$ DDM", zorder=1)
ax_Y.plot(strainax(F), Eij_dsc[:,I11], '--', c='tab:green', label=r"$E_{11}$ DDM", zorder=1)

if SHOW_ALT_EIJ_COMPONENTS:
    ax_Y.plot(strainax(F), Eij_dsc[:,I22], '-.', c='tab:purple', label=r"$E_{22}$ DDM", zorder=1)
    ax_Y.plot(strainax(F), Eij_dsc[:,I33], ':',  c='tab:pink',   label=r"$E_{33}$ DDM", zorder=1)

# Axis labels

F_indices = ('zz' if T_EXP==T_EXP_UC else 'xz')
xaxlbl = r'\gamma_{%s}'%(F_indices) if T_EXP==T_EXP_SS else r'\epsilon_{%s}'%(F_indices)
ax_Y.set_xlabel('$%s$'%(xaxlbl))
ax_Y.set_ylabel(r'$E_{ij}$')

# Legend

legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':1, 'ncol':1, 'handlelength':1.34, 'labelspacing':0.25, 'fontsize':FSSMALL}
leg = ax_Y.legend(loc=1, bbox_to_anchor=((1,1) if T_EXP==T_EXP_UC else (1,0.9)), **legkwargs)
leg.get_frame().set_linewidth(0.8);
sfplt.panellabel(ax_Y, 2, r'\textit{(a)}', frameon=False, fontsize=FS+2, bbox=(-0.15,1.15))

# Axis ticks and limits

if T_EXP==T_EXP_SS: xticks = np.arange(0,20,0.25) 
if T_EXP==T_EXP_UC: xticks = np.arange(0,1,0.1)
ax_Y.set_xticks(xticks, minor=True)  
ax_Y.set_xlim(xlims)

ax_Y.set_yticks(np.arange(0,10+1e-3,0.5))
ax_Y.set_yticks(np.arange(0,10+1e-3,0.1), minor=True)
ax_Y.set_ylim(ylims)

### Plot ODFs

for ii,nn in enumerate(ODF_tsteps):

    ax_blm = ax_ODF[ii*2]
    ax_blm.set_global()
    
    ax_nlm = ax_ODF[ii*2+1]
    ax_nlm.set_global()
    
    sfplt.plotODF(blm[nn,:], lm, ax_blm, cblabel=r'$b/N$', cmap='Reds',  lvlset=[np.linspace(0,0.41,9), lambda x,p:'%.1f'%x], nchunk=None)
    sfplt.plotODF(nlm[nn,:], lm, ax_nlm, cblabel=r'$n/N$', cmap='Blues', lvlset=[np.linspace(0,0.41,9), lambda x,p:'%.1f'%x], nchunk=None)

    ax_blm.set_title(r'$%s=%.1f$'%(xaxlbl, strainax(F)[nn]), fontsize=FS)
    ax_nlm.set_title(r'$%s=%.1f$'%(xaxlbl, strainax(F)[nn]), fontsize=FS)

    sfplt.panellabel(ax_blm, 2, r'\textit{(%s)}'%(chr(ord('b')+ii)), frameon=False, fontsize=FS+2, bbox=(-0.75,1.6))

    kwargs = dict(transform=geo, fontsize=FSSMALL, color='k')
    for ii, ei in enumerate((e1[nn,:],e2[nn,:],e3[nn,:]),1):
        text = r'$\vb{m}_{%i}$'%(ii)
        sfplt.plotS2text(ax_blm, +ei, text, **kwargs)
        sfplt.plotS2text(ax_blm, -ei, text, **kwargs)
        sfplt.plotS2text(ax_nlm, +ei, text, **kwargs)
        sfplt.plotS2text(ax_nlm, -ei, text, **kwargs)

### Plot parcel deformation as inset

if T_EXP==T_EXP_SS: parcel_tsteps = ODF_tsteps[[0,]]
if T_EXP==T_EXP_UC: parcel_tsteps = ODF_tsteps[[1,]]

for ii,nn in enumerate(parcel_tsteps):

    # Parcel locations (normalized figure coords, so must be set manually like this)
    if T_EXP==T_EXP_SS: 
        y0 = 0.76
        pc = np.array([[0.21,y0],])
        
    if T_EXP==T_EXP_UC: 
        y0 = 0.77
        pc = np.array([[0.475,y0]])

    axs = 0.15
    ax_sub = fig.add_axes([pc[ii,0], pc[ii,1], axs, axs], projection='3d')
    ax_sub.patch.set_alpha(0.0)
    lw = 0.7
    sfplt.plotparcel(ax_sub, F[nn,:,:], lw=lw, lwinit=lw, fonttex=True)
    ax_sub.set_title(r'$%s=%.1f$'%(xaxlbl,F[nn,i,j]), fontsize=FSSMALL,  x=0.5, y=0.92)
        
### Save figure

T_EXP_STR = {T_EXP_UC:'UC', T_EXP_SS:'SS'}
fname = 'Eij-lin-sachs-SDM-vs-DDM--%s.pdf'%(T_EXP_STR[T_EXP])
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=175)

