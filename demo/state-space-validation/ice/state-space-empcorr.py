# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from scipy.optimize import curve_fit
import pickle, glob

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

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

L=8
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
# Determine ideal boundary line
#--------------------

m = [0,0,1]
Il24 = [sf.I20, sf.I40] # l=2,4, m=0 coefs

# delta distributed
n20_unidir, n40_unidir = np.real(sf.nlm_ideal(m, 0, L))[Il24]/normfac

# x--y planar distributed
n20_planar, n40_planar = np.real(sf.nlm_ideal(m, np.pi/2, L))[Il24]/normfac    

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

#n20_delta, n20_planar = np.sqrt(5), -np.sqrt(5)/2 # 2.23606797749979, -1.1180340067479093
print(n20_planar,n20_unidir)
x_corr = np.linspace(n20_planar,n20_unidir,100)
p_model = np.poly1d(np.flipud(nhat40_pc))
#y_corr = p_model(x_corr) # debug
y_corr = sf.nhat40_empcorr_ice(x_corr) 

#--------------------
# IBOF closure from Elmer/Ice
#--------------------

y_IBOF = y_corr * 0 # init

for ii, x in enumerate(x_corr):
    nlm_hat = np.zeros((nlm_len), dtype=np.complex128)
    nlm_hat[[0,sf.I20]] = [1, x]
    nlm_ = sf.a4_to_nlm(sf.a4_IBOF(sf.a2(nlm_hat*1/np.sqrt(4*np.pi))))
    y_IBOF[ii] = np.real(nlm_[sf.I40]/nlm_[0]) # nhat40 

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

print('Determining subspace of valid eigenvalues...')
xv, yv = np.meshgrid(x, y, indexing='xy')
validregion = np.reshape(sf.nlm_isvalid(xv.flatten(), yv.flatten()), (RESY, RESX))

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

# Unidirectional
ax.plot(n20_unidir,n40_unidir, marker='o', ms=mse, ls='none', c=c_smax, label=None, zorder=20)
plt.text(n20_unidir-0.1, n40_unidir+dytext, '{\\bf Unidirectional}', color=c_smax, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# Planar
ax.plot(n20_planar, n40_planar, marker='o', ms=mse, ls='none', c=c_girdle, label=None, zorder=20)
plt.text(n20_planar, n40_planar+dytext, '{\\bf Planar}', color=c_girdle, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# Shading labels
plt.text(0.26/normfac, 0.1/normfac, '{\\bf Single maximum}', color=c_smax, ha='center', rotation=30, fontsize=FSANNO)
plt.text(-0.175/normfac, 0.125/normfac, '{\\bf Girdle}', color=c_girdle, ha='center', rotation=-40, fontsize=FSANNO)
plt.text(2.1, -1, '{\\bf Unphysical}\n\\bf{eigenvalues}', color='0.3', ha='center', rotation=0, fontsize=FSANNO)

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

ax.plot(x_corr, y_IBOF, ':k', lw=2, label='IBOF closure', zorder=10)

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
#secax.set_xlabel(r'$a^{(2)}_{zz}$')
secax.set_xlabel(r'$\ev*{c_i c_j}_{zz}$')
xticks = np.arange(0,1+1e-3,0.1)
secax.set_xticks(xticks[0::2])
secax.set_xticks(xticks[::1], minor=True)


### Plot ODF insets?

if 1:

    geo, prj = sfplt.getprojection(rotation=55+180, inclination=50)

    W = 0.125 # ax width
    
    ### Modeled data
    
    arrmag = 0.125
    arr = lambda ang: arrmag*np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
    n00 = 1/np.sqrt(4*np.pi)
    nlm_constructor = lambda nhat20: np.array([n00, 0,0,nhat20*n00,0,0, 0,0,0,0,sf.nhat40_empcorr_ice(nhat20)[0]*n00,0,0,0,0], dtype=np.complex128)
    ODF_plots = (\
        {'nlm':nlm_constructor(-0.75), 'title':'', 'axloc':(0.205, 0.16), 'darr':arr(-90), 'lvlmax':0.4}, \
        {'nlm':nlm_constructor(+1.00), 'title':'', 'axloc':(0.577, 0.17), 'darr':arr(-90), 'lvlmax':0.4}, \
        {'nlm':nlm_constructor(+1.75), 'title':'', 'axloc':(0.737, 0.33), 'darr':arr(-90), 'lvlmax':0.4}, \
    )

    for ODF in ODF_plots:

        axpos = [ODF['axloc'][0],ODF['axloc'][1], W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        nlm = ODF['nlm']
        lvlset = [np.linspace(0.0,ODF['lvlmax'],6), lambda x,p:'%.1f'%x]
        sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, showcb=False)

        # Arrow to ODF state
        n20_, n40_ = np.real(nlm[sf.I20])/normfac, np.real(nlm[sf.I40])/normfac
        ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                        xytext=(n20_+ODF['darr'][0]/normfac, n40_+sc**2*ODF['darr'][1]/normfac), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", linewidth=1.5, edgecolor='0.2', facecolor='0.2'),zorder=20)            

        axin.set_title(ODF['title'], fontsize=FS-3)

### Save figure

#plt.tight_layout()
plt.savefig('state-space-empcorr.png', transparent=1,  dpi=175)
plt.close()

