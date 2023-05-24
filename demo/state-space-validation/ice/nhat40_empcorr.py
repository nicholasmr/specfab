# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from scipy.optimize import curve_fit
import pickle, glob

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)
sys.path.insert(0, '../../../demo') # for importing local specfabpy build (if available) and common python header
from header import * # contains matplotlib setup etc.
from specfabpy import specfabpy as sf

SELFNAME = 'state-space-validation' # used as prefix for pickled files
os.system('mkdir -p specfab-state-trajectories')
def pfile(fname): return "specfab-state-trajectories/%s--%s.p"%(SELFNAME, fname) # full path name for pickled files

#--------------------
# Flags
#--------------------

DEBUG = 0 # Low-resolution plotting, etc.
            
#--------------------
# Config
#--------------------

if DEBUG: 
    RES = 2*50 
    RESX, RESY = int(RES/2), RES
else:
    RESX = RESY = 5*100
    
experiments = (expr_GRIP, expr_LAWDOME,  expr_EGRIP_MAIN, expr_SPICE, expr_Priestley, expr_EGRIP_S5) 

L=4
lm, nlm_len = sf.init(L) 

#--------------------
# Load pre-processed correlations (see process-experiments.py)
#--------------------

correlations = []

for expr in experiments:
    fcorr = "observed-states/%s.p"%(expr['path']) 
    print('=== Loading correlations: %s ==='%(fcorr))
    corr = pickle.load(open(fcorr, "rb"))
    corr['n20'] /= normfac # Normalization
    corr['n40'] /= normfac
    corr['n60'] /= normfac
    correlations.append(corr) 

#--------------------
# Polynomial-fitted correlation curve used by Gerber et al. (2022)
#--------------------

x_, y_ = [], []
for ii,expr in enumerate(experiments):  
    x_.append(correlations[ii]['n20']) # n_2^0
    y_.append(correlations[ii]['n40']) # n_4^0
x = np.concatenate(x_)
y = np.concatenate(y_)

#def a function
def f(x, a, b, c, d):
    return a*x + b*x**2 + c*x**3 + d*x**4 #+ e*x**5

nhat40_pc, pcov = curve_fit(f, x, y)
nhat40_pc = np.concatenate(([0], nhat40_pc))
print("Correlation coefs in f(x) = a*x + b*x**2 + c*x**3 + d*x**4 are\n", nhat40_pc)

n20_delta, n20_planar = np.sqrt(5), -np.sqrt(5)/2 # 2.23606797749979, -1.1180340067479093
x_corr = np.linspace(n20_planar,n20_delta,100)
p_model = np.poly1d(np.flipud(nhat40_pc))
#y_corr = p_model(x_corr) # debug
y_corr = sf.nhat40_empcorr_ice(x_corr) 

#--------------------
# Construct plot
#--------------------

ms = 8
mse = ms+0.5 # end-member case points
FSLEG = FS-0.5
FSANNO = FS-1.5

c_girdle = '#8c510a' 
c_smax   = '#01665e'
c_cdrx   = 'k'

cl_smax   = np.array([199,234,229])/255
cl_girdle = np.array([246,232,195])/255

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'columnspacing': 0.5, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

### Setup figure

scale = 4.0
fig = plt.figure(figsize=(1.3*scale,1.1*scale))
plt.subplots_adjust(left=0.12, right=0.99, top=0.88, bottom=0.13)
ax = plt.gca()
xlims, ylims = [-1.45,2.65], [-1.4,3.5]
sc = np.diff(ylims)/np.diff(xlims)
x = np.linspace(xlims[0],xlims[1],RESX)
y = np.linspace(ylims[0],ylims[1],RESY)

### Determine valid subspace (valid eigenvalues)

validregion = np.zeros((RESY, RESX)) # 0 = invalid, 1 = valid 
validregion_lowerbound = np.zeros((RESX)) # points along lower boundary, used for colouring the background (shading) of subspace with ~circle fabrics.
print('Determining subspace of valid eigenvalues...', end='')
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        nlm_ = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients
        nlm_[0], nlm_[3], nlm_[10] = 1, x_, y_
        a2_ = sf.a2(nlm_) # diagional
        a2_eigvals = np.sort(np.diag(a2_))
        Q1,Q2,Q3,Q4,Q5,Q6, a4_eigvals = sf.a4_eigentensors(nlm_)
        isvalid_a2_eig = (np.amin(a2_eigvals)>=0) and (np.amax(a2_eigvals)<=1)  
        isvalid_a4_eig = (np.amin(a4_eigvals)>=0) and (np.amax(a4_eigvals)<=1)
        validregion[yii,xii] = isvalid_a2_eig and isvalid_a4_eig
    validregion_lowerbound[xii] = y[np.argmax(validregion[:,xii])]
            
print('done')

### Determine subspace shadings

imdat = np.empty((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)

print('Determining shading regions for different fabric types...', end='')
I = np.nonzero(np.abs(x_corr) > 0.4)
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        if validregion[yii,xii]:
            distnn = np.amin( np.sqrt(np.real((x_-x_corr[I])**2 + (1/sc*(y_-y_corr[I]))**2)) ) # distance from model-line points
            var, expo = 1e-3, 6
            imdat[yii,xii,-1] = np.exp(-distnn**expo/var) # set alpha depending on distance to model lines/points
            if x_<0: imdat[yii,xii,0:-1] = cl_girdle
            else:    imdat[yii,xii,0:-1] = cl_smax
        else: 
            imdat[yii,xii,0:-1] = [0.85]*3 # bad 
            imdat[yii,xii,-1] = 1
print('done')

# Plot valid/invalid subspaces and fabric-type shadings within valid subspace
imdat0 = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
im = ax.imshow(imdat0, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
im = ax.imshow(imdat, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)

### End-member cases

# Isotropic state
ax.plot(0,0,'o', ms=mse, c='k', label=None, zorder=20)
dytext = 0.04/normfac
plt.text(0, 0+dytext, r'{\bf Isotropic}', color='k', ha='center', va='bottom', fontsize=FSANNO)

# Unidirectional/delta-function (single max)
n20_delta = np.real(sp.sph_harm(0, 2, 0,0))/normfac
n40_delta = np.real(sp.sph_harm(0, 4, 0,0))/normfac
ax.plot(n20_delta,n40_delta, marker='o', ms=mse, ls='none', c=c_smax, label=None, zorder=20)
plt.text(n20_delta-0.1, n40_delta+dytext, '{\\bf Unidirectional}', color=c_smax, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# Planar (great circl) isotropy 
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
ax.plot(x_, y_, marker='o', ms=mse, ls='none', c=c_girdle, label=None, zorder=20)
plt.text(x_, y_+dytext, '{\\bf Planar}\n\\bf{isotropic}', color=c_girdle, ha='center', va='bottom', ma='center', fontsize=FSANNO)

print(n40_delta, y_)

# Shading labels
plt.text(0.26/normfac, 0.1/normfac, '{\\bf Single maximum}', color=c_smax, ha='center', rotation=30, fontsize=FSANNO)
plt.text(-0.175/normfac, 0.125/normfac, '{\\bf Girdle}', color=c_girdle, ha='center', rotation=-40, fontsize=FSANNO)

plt.text(-0.25/normfac, -0.32/normfac, '{\\bf Unphysical}\n\\bf{eigenvalues}', color='0.3', ha='center', rotation=0, fontsize=FSANNO)

### Experimental data points

for ii,expr in enumerate(experiments):
    x = correlations[ii]['n20'] # n_2^0
    y = correlations[ii]['n40'] # n_4^0
    if expr['type']=='ss':  mrk = 's'
    if expr['type']=='ue':  mrk = '^'
    if expr['type']=='uc':  mrk = 'd'
    if expr['type']=='ucw': mrk = 'X'
    ax.plot(x,y, ms=ms, ls='none', color=expr['color'], fillstyle='none', marker=mrk, label=expr['plotname'], zorder=1+(len(experiments)-ii))

### Plot correlation

ax.plot(x_corr, y_corr, '-k', lw=2, label='Empirical correlation', zorder=10)

### Aux

plt.sca(ax)
plt.xlabel(r'$\hat{n}_2^0$')
plt.ylabel(r'$\hat{n}_4^0$')

leg = plt.legend(loc=2, fontsize=FSLEG, frameon=False, ncol=2, **legkwargs); 

### Limits

plt.xlim(xlims)
plt.ylim(ylims)

### Second x axis for a_zz^(2) comparrison 

secax = ax.secondary_xaxis('top', functions=(n20_to_azz, azz_to_n20))
secax.set_xlabel(r'$a^{(2)}_{zz}$')
xticks = np.arange(0,1+1e-3,0.1)
secax.set_xticks(xticks[0::2])
secax.set_xticks(xticks[::1], minor=True)

### Save figure

#plt.tight_layout()
plt.savefig('nhat40_empcorr_ice.png', transparent=1,  dpi=175)
plt.close()

