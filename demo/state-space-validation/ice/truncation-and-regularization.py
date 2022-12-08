# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Fabric state-space diagram quantifying the effect of truncation (L) and (hyper)regularization (Lilien et al., 2022)
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import pickle
from progress.bar import Bar

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from localheader import *
sys.path.insert(0, '../../../demo') # for importing local specfabpy build (if available) and common python header
from header import * # contains matplotlib setup etc.
from specfabpy import specfabpy as sf

# This is the presumed psi_0^0 coefficient when deriving nlm from a^(2) or a^(4), a central normalization factor used below. 
normfac = 1/np.sqrt(4*np.pi) 

SELFNAME = sys.argv[0][:-3] # used as prefix for pickled files
os.system('mkdir -p specfab-state-trajectories')
def pfile(fname): return "specfab-state-trajectories/%s--%s.p"%(SELFNAME, fname) # full path name for pickled files

#--------------------
# Flags
#--------------------

INTEGRATE_MODEL = 1  # Generate model lines from scratch? Else load saved Pickle files.
DEBUG           = 0  # For faster plotting (lower resolution)
            
#--------------------
# Config
#--------------------

Nt = 300 # Number of integration steps for model trajectories (should be large enough to reach steady state)
RESX = RESY = 1*100 if DEBUG else 15*100 # Fast plotting (debug)?

#--------------------
# Modeled correlations
#--------------------

def integrate_model(nlm0, ugrad, dt, reg='hyper', Nt=Nt, rotate=False, name=None):

    eps = (ugrad+np.transpose(ugrad))/2 # Symmetric part (strain-rate)
    omg = (ugrad-np.transpose(ugrad))/2 # Anti-symmetric part (spin)

    nlm = np.zeros((Nt,nlm_len), dtype=np.complex64)
    nlm[0,:] = nlm0 # Initial state 
             
    # Euler integration 
    with Bar('dt=%.3e, Nt=%i :: L=%i (nlm_len=%i) ::'%(dt,Nt,L,nlm_len), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds') as bar:
        for tt in np.arange(1,Nt):
            nlm_prev = nlm[tt-1,:]
            iota, zeta = 1,0
            M_LROT = sf.M_LROT(nlm_prev, eps, omg, iota, zeta)
            M_REG  = sf.M_REG(nlm_prev, eps) if reg=='hyper' else 10e-2*sf.M_CDRX(nlm_prev) # note that CDRX is identical the Laplacian operator, so used as regularization here
            M      = M_LROT + M_REG
            nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev)
            nlm[tt,:] = sf.apply_bounds(nlm[tt,:])
            bar.next()

    if name is not None: pickle.dump([nlm, lm, nlm_len], open(pfile(name), "wb"))

    return nlm


if INTEGRATE_MODEL:

    print('*** Generating model trajectories from scratch. You can re-use the saves trajectories to avoid re-calculating them by setting INTEGRATE_MODEL=1')

    ### Solve for state-vector time evolution
    
    L = 8 if DEBUG else 20 # Example of high-order truncation
    lm, nlm_len = sf.init(L) 
    nlm_iso = np.zeros((nlm_len), dtype=np.complex64)
    nlm_iso[0] = normfac
    nlm_uc20 = integrate_model(nlm_iso, +np.diag([.5, .5, -1]), 0.015*1, name='uc20') # unconfined compression
    nlm_ue20 = integrate_model(nlm_iso, -np.diag([.5, .5, -1]), 0.007*2, name='ue20') # unconfined extension 
    
    L=6 # Example of low-order truncation
    lm, nlm_len = sf.init(L) 
    nlm_iso = np.zeros((nlm_len), dtype=np.complex64)
    nlm_iso[0] = normfac
    nlm_uc6 = integrate_model(nlm_iso, +np.diag([.5, .5, -1]), 0.015*1, name='uc6') # unconfined compression 
    nlm_ue6 = integrate_model(nlm_iso, -np.diag([.5, .5, -1]), 0.007*2, name='ue6') # unconfined extension 

    nlm_uc6_lap = integrate_model(nlm_iso, +np.diag([.5, .5, -1]), 0.015*1, name='uc6lap', reg='lap') # unconfined compression 
    nlm_ue6_lap = integrate_model(nlm_iso, -np.diag([.5, .5, -1]), 0.007*2, name='ue6lap', reg='lap') # unconfined extension 


### Load solutions
nlm_uc20,   lm, nlm_len = pickle.load(open(pfile('uc20'), "rb"))
nlm_ue20,   lm, nlm_len = pickle.load(open(pfile('ue20'), "rb"))
nlm_uc6,    lm, nlm_len = pickle.load(open(pfile('uc6'), "rb"))
nlm_ue6,    lm, nlm_len = pickle.load(open(pfile('ue6'), "rb"))       
nlm_uc6lap, lm, nlm_len = pickle.load(open(pfile('uc6lap'), "rb"))
nlm_ue6lap, lm, nlm_len = pickle.load(open(pfile('ue6lap'), "rb"))       

### Normalize
nlm_uc20   = np.array([ nlm_uc20[tt,:]/nlm_uc20[tt,0] for tt in np.arange(Nt) ])
nlm_ue20   = np.array([ nlm_ue20[tt,:]/nlm_ue20[tt,0] for tt in np.arange(Nt) ])
nlm_uc6    = np.array([ nlm_uc6[tt,:]/nlm_uc6[tt,0]   for tt in np.arange(Nt) ])
nlm_ue6    = np.array([ nlm_ue6[tt,:]/nlm_ue6[tt,0]   for tt in np.arange(Nt) ])
nlm_uc6lap = np.array([ nlm_uc6lap[tt,:]/nlm_uc6lap[tt,0] for tt in np.arange(Nt) ])
nlm_ue6lap = np.array([ nlm_ue6lap[tt,:]/nlm_ue6lap[tt,0] for tt in np.arange(Nt) ])

#--------------------
# Construct plot
#--------------------

ms = 6.0 
mse = 7 # point size of end-member cases
FSLEG = FS-1.5
FSANNO = FS-2.5

c_girdle = '#8c510a' 
c_smax   = '#01665e'
c_girdlelight = '#dfc27d'
c_smaxlight   = '#7fcdbb'

legkwargs = {'handlelength':1.4, 'framealpha':0.7, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

### Setup figure

scale = 3.5
fig = plt.figure(figsize=(1*1.55*scale,1.5*scale))
gs = gridspec.GridSpec(1, 1, hspace=0)
plt.subplots_adjust(left=0.12, right=0.98, top=0.9, bottom=0.3)
ax_Exz  = fig.add_subplot(gs[0,0])
xlims, ylims = [-1.45,2.65], [-1.8,3.75]
sc = np.diff(ylims)/np.diff(xlims)
x = np.linspace(xlims[0],xlims[1],RESX)
y = np.linspace(ylims[0],ylims[1],RESY)
X,Y = x.copy(), y.copy()

### Determine E_zz enhancement factor

Ezz = np.zeros((RESY, RESX)) *np.nan
Exz = np.zeros((RESY, RESX)) *np.nan
Epq = np.zeros((RESY, RESX)) *np.nan

# Monocrystal fluid parameters for enhancement-factor calculation
Ecc,nprime = 1, 1
Eca,alpha = 1e3, 0.0125

# Symmetry axis "m" and transverse axis "t"
m, t = np.array([0,0,1]), np.array([1,0,0])
p, q = (m+t)/np.sqrt(2), (m-t)/np.sqrt(2)
mm, mt, pq = np.tensordot(m,m, axes=0), np.tensordot(m,t, axes=0), np.tensordot(p,q, axes=0)

# Idealized stress states (aligned with fabric axes)
tau_ps_mm = 1*(np.identity(3)-3*mm) 
tau_ss_mt = 1*(mt + np.transpose(mt)) 
tau_ss_pq = 1*(pq + np.transpose(pq))

print('Calculating E_zz map ...', end='')
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        nlm_ = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients
        nlm_[0], nlm_[3], nlm_[10] = 1, x_, y_
        Ezz[yii,xii] = sf.Evw(nlm_, mm,tau_ps_mm, Ecc,Eca,alpha,nprime)
        Exz[yii,xii] = sf.Evw(nlm_, mt,tau_ss_mt, Ecc,Eca,alpha,nprime)

print('done')

# Verify that end-member case of perfectly aligned c-axes (single max or girdle) reproduces monocrystal behavior
if DEBUG:
    a2 = np.tensordot(m,m, axes=0)
    a4 = np.tensordot(a2,a2, axes=0)
    nlm_sm = sf.a4_to_nlm(a4)
    Ezz_sm = sf.Evw(nlm_sm, mm,tau_ps_mm, Ecc,Eca,alpha,nprime)
    Exz_sm = sf.Evw(nlm_sm, mt,tau_ss_mt, Ecc,Eca,alpha,nprime)
    Epq_sm = sf.Evw(nlm_sm, pq,tau_ss_pq, Ecc,Eca,alpha,nprime)
    print(nlm_sm)
    print('Ezz_sm, Exz_sm, Exz_sm/Epq_sm = %.1e, %.1e, %.1e'%(Ezz_sm,Exz_sm, Exz_sm/Epq_sm))

### More plot setup

dlvl = 1
lvls = np.arange(0,6,dlvl)
cmap = plt.cm.Blues.copy()
cmap.set_under('#d6604d')
CS = ax_Exz.contourf(x,y,Exz, levels=lvls, cmap=cmap, extend='both')
hcb = plt.colorbar(CS, ax=ax_Exz, orientation='vertical') # , aspect=25
hcb.set_label(r'$E_{xz}$')

lvls = [0.01, 0.1, 0.5, 1, 1.5, 2]
CS = ax_Exz.contour(x,y,Ezz, levels=lvls, linewidths=0.7, linestyles='-', colors='k')
manual_locations = [(-0.2, 1.5)] + [(0.51, y_) for y_ in np.linspace(-1.4, 2.1, len(lvls[1:]))]
ax_Exz.clabel(CS, CS.levels, inline=True, fmt=r'$E_{zz}=%.2f$', fontsize=FS-1, manual=manual_locations)

ax = ax_Exz

### Model lines

h_uc20   = plot_trajectory(ax, nlm_uc20,   arrpos=9,    c=c_smax,        ls='--', endmarker=True, mse=mse)
h_ue20   = plot_trajectory(ax, nlm_ue20,   arrpos=17,   c=c_girdle,      ls='--', endmarker=True, mse=mse)
h_uc6    = plot_trajectory(ax, nlm_uc6,    arrpos=None, c=c_smax,        ls='-',  endmarker=True, mse=mse)
h_ue6    = plot_trajectory(ax, nlm_ue6,    arrpos=None, c=c_girdle,      ls='-',  endmarker=True, mse=mse)
h_uc6lap = plot_trajectory(ax, nlm_uc6lap, arrpos=None, c=c_smaxlight,   ls=':',  endmarker=True, mse=mse)
h_ue6lap = plot_trajectory(ax, nlm_ue6lap, arrpos=None, c=c_girdlelight, ls=':',  endmarker=True, mse=mse)

h_modellines = [h_ue20, h_ue6, h_ue6lap,   h_uc20, h_uc6, h_uc6lap]
legend_strings = ['Unconf. compr., L=20', 'Unconf. compr., L=6', 'Unconf. compr., L=6, no hyper dif.', \
                  'Unconf. exten., L=20', 'Unconf. exten., L=6', 'Unconf. exten., L=6, no hyper dif.']

legend_modellines = ax.legend(h_modellines, legend_strings, title=r'{\bf Modelled fabric state trajectories}', title_fontsize=FSLEG, \
    loc=3, bbox_to_anchor=(-0.17,-0.5), ncol=2, fontsize=FS-1, frameon=True, **legkwargs)

### End-member cases

# Isotropic state
ax.plot(0,0,'o', ms=mse, c='k', label=None, zorder=20)
dytext = 0.04/normfac
ax.text(0, 0-dytext, r'{\bf Isotropic}', color='k', ha='center', va='top', fontsize=FSANNO)

# Unidirectional/delta-function (single max)
n20_delta = np.real(sp.sph_harm(0, 2, 0,0))/normfac
n40_delta = np.real(sp.sph_harm(0, 4, 0,0))/normfac
#print(n20_delta*normfac,n40_delta*normfac)
ax.plot(n20_delta,n40_delta, marker='o', ms=mse, ls='none', c=c_smax, label=None)
ax.text(n20_delta-0.4*dytext, n40_delta+dytext, '{\\bf Unidirectional}', color=c_smax, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# Planar (great circle) isotropy 
x, y = np.array([1,0,0]), np.array([0,1,0])
x2, y2 = np.einsum('i,j',x,x), np.einsum('i,j',y,y)
x4, y4 = np.einsum('i,j,k,l',x,x,x,x), np.einsum('i,j,k,l',y,y,y,y)
xy_sym = np.einsum('i,j',x,y) + np.einsum('i,j',y,x)
xy2 = x2 + y2
a2 = x2/2 + y2/2
a4 = x4/4 + y4/4 + np.einsum('ij,kl',xy2,xy2)/8 + np.einsum('ij,kl',xy_sym,xy_sym)/8 
nlm_girdle = np.real(sf.a4_to_nlm(a4))
#print(nlm_girdle, nlm_girdle[3],nlm_girdle[10])
x_, y_ = np.real(nlm_girdle[3])/normfac, np.real(nlm_girdle[10])/normfac
ax.plot(x_, y_, marker='o', ms=mse, ls='none', c=c_girdle, label=None)
ax.text(x_, y_+dytext, '{\\bf Planar}\n\\bf{isotropic}', color=c_girdle, ha='center', va='bottom', ma='center', fontsize=FSANNO)

### Power-spectrum cap
if 1:
    badspec = np.ones((RESY, RESX)) 
    I = np.argmin(np.abs(X-n20_delta))
    badspec[:,I:] = 0
    I = np.argmin(np.abs(Y-n40_delta))
    badspec[I:,:] = 0
    lvls = [0, 0.5]
    plt.rcParams['hatch.color'] = '0.3'
    plt.rcParams['hatch.linewidth'] = 0.6
    cs = ax.contourf(X, Y, badspec, hatches=["\\",None], levels=lvls, cmap='gray', alpha=0)
    cs = ax.contour( X, Y, badspec, levels=[0.5,], linewidths=plt.rcParams['hatch.linewidth'], colors=plt.rcParams['hatch.color'])
    ax.text(-1.2, 3.3, r'\bf Restricted', fontsize=FS-2, color='0.1')
else:
    lw = 1.1
    ax.plot([n20_delta,n20_delta],[-10,n40_delta],':k', lw=lw)
    ax.plot([-10,n20_delta],[n40_delta,n40_delta],':k', lw=lw)

### Aux
ax.set_xlabel(r'$\hat{\psi}_2^0$')
ax.set_ylabel(r'$\hat{\psi}_4^0$')

### Limits
ax.set_xlim(xlims)
ax.set_ylim(ylims)

### Second x axis
secax = ax.secondary_xaxis('top', functions=(n20_to_azz, azz_to_n20))
secax.set_xlabel(r'$a^{(2)}_{zz}$')

### Print quantification of truncation error
if DEBUG:
    n20_delta = np.real(sp.sph_harm(0, 2, 0,0))/normfac
    n40_delta = np.real(sp.sph_harm(0, 4, 0,0))/normfac
    print(nlm_uc20[-1,3]/n20_delta, nlm_uc20[-1,10]/n40_delta)
    print( nlm_uc6[-1,3]/n20_delta,  nlm_uc6[-1,10]/n40_delta)

### Save figure
plt.savefig('%s.png'%(SELFNAME), dpi=175)
plt.close()

