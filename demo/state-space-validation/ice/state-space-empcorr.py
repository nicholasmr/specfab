# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from scipy.optimize import curve_fit
import pickle, glob

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors

from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)

sys.path.append('../../')
import demolib as dl

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-0.5
FSANNO = FS-1.5

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

L = 8
lm, nlm_len = sf.init(L) 

#--------------------
# Load pre-processed correlations (see process-experiments.py)
#--------------------

correlations = []

for expr in experiments:
    fcorr = "observed-states/%s.p"%(expr['path']) 
    print('=== Loading correlations: %s ==='%(fcorr))
    corr = pickle.load(open(fcorr, "rb"))
    corr['n20'] /= norm # Normalization
    corr['n40'] /= norm
    corr['n60'] /= norm
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

def f(x, a, b, c, d):
    return a*x + b*x**2 + c*x**3 + d*x**4 #+ e*x**5

nhat40_pc, pcov = curve_fit(f, x, y)
nhat40_pc = np.concatenate(([0], nhat40_pc))
print("Correlation coefs in f(x) = a*x + b*x**2 + c*x**3 + d*x**4 are\n", nhat40_pc)

nl0_unidir, nl0_planar, nl0_circle = sfdsc.nlm_ideal_cases(norm=norm)
x_corr = np.linspace(nl0_planar[0],nl0_unidir[0],100)
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
    nlm_ = sf.a4_to_nlm(sf.a4_IBOF(sf.a2(nlm_hat*norm)))
    y_IBOF[ii] = np.real(nlm_[sf.I40]/nlm_[0]) # nhat40 

#--------------------
# Construct plot
#--------------------

### Setup figure

scale = 4.0
fig = plt.figure(figsize=(1.3*scale,1.1*scale))
plt.subplots_adjust(left=0.12, right=0.99, top=0.88, bottom=0.13)
ax = plt.gca()

xlims, ylims = [-1.45,2.65], [-1.4,3.5]
sc = np.diff(ylims)/np.diff(xlims)

### Determine valid subspace (valid eigenvalues)

print('Determining subspace of valid eigenvalues...')

isvalid, x,y = dl.nlm_isvalid_grid(xlims, ylims, RESX, RESY)

print('Determining shading regions for different fabric types...', end='')
imdat = np.empty((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
I = np.nonzero(np.abs(x_corr) > 0.4)
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        if isvalid[yii,xii]:
            distnn = np.amin( np.sqrt(np.real((x_-x_corr[I])**2 + (1/sc*(y_-y_corr[I]))**2)) ) # distance from model-line points
            var, expo = 1e-3, 6
            imdat[yii,xii,-1] = np.exp(-distnn**expo/var) # set alpha depending on distance to model lines/points
            if x_<0: imdat[yii,xii,0:-1] = matplotlib.colors.to_rgb(dl.cvl_planar)
            else:    imdat[yii,xii,0:-1] = matplotlib.colors.to_rgb(dl.cvl_unidir)
        else: 
            imdat[yii,xii,:-1] = matplotlib.colors.to_rgb(sfplt.c_lgray) # bad 
            imdat[yii,xii,-1] = 1
print('done')

# Plot valid/invalid subspaces and fabric-type shadings within valid subspace
imdat0 = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
kwargs = dict(aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
im = ax.imshow(imdat0, **kwargs)
im = ax.imshow(imdat,  **kwargs)

### End-member cases

dl.plot_nlm_cases(ax, FSANNO, ms=8.5, dy0=0.07, show_circle=False)

# Shading labels
plt.text(0.26/norm, 0.1/norm, '{\\bf Single maximum}', color=dl.c_unidir, ha='center', rotation=30, fontsize=FSANNO)
plt.text(-0.175/norm, 0.125/norm, '{\\bf Girdle}', color=dl.c_planar, ha='center', rotation=-40, fontsize=FSANNO)
plt.text(2.1, -1, '{\\bf Unphysical}\n\\bf{eigenvalues}', color=sfplt.c_dgray, ha='center', rotation=0, fontsize=FSANNO)

### Experimental data points

for ii,expr in enumerate(experiments):
    x = correlations[ii]['n20'] # n_2^0
    y = correlations[ii]['n40'] # n_4^0
    if expr['type']=='ss':  mrk = 's'
    if expr['type']=='ue':  mrk = '^'
    if expr['type']=='uc':  mrk = 'd'
    if expr['type']=='ucw': mrk = 'X'
    ax.plot(x,y, ms=8, ls='none', color=expr['color'], fillstyle='none', marker=mrk, label=expr['plotname'], zorder=1+(len(experiments)-ii))

### Plot correlation lines

ax.plot(x_corr, y_corr, '-k', lw=2, label='Empirical correlation', zorder=10)
ax.plot(x_corr, y_IBOF, ':k', lw=2, label='IBOF closure', zorder=10)

### Misc

plt.sca(ax)

plt.xlabel(r'$\hat{n}_2^0$')
plt.ylabel(r'$\hat{n}_4^0$')

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'columnspacing': 0.5, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}
leg = plt.legend(loc=2, fontsize=FSLEG, frameon=False, ncol=2, **legkwargs); 

plt.xlim(xlims)
plt.ylim(ylims)

### Second x axis for a_zz^(2) comparrison 

secax = ax.secondary_xaxis('top', functions=(n20_to_azz, azz_to_n20))
secax.set_xlabel(r'$\ev*{c_i c_j}_{zz}$')
xticks = np.arange(0,1+1e-3,0.1)
secax.set_xticks(xticks[0::2])
secax.set_xticks(xticks[::1], minor=True)

### Plot ODF insets?

if 1:

    geo, prj = sfplt.getprojection(rotation=55+180, inclination=50)

    ### Modeled data
    
    arr = lambda ang: 0.125*np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
    n00 = norm
    nlm_constructor = lambda nhat20: np.array([n00, 0,0,nhat20*n00,0,0, 0,0,0,0,sf.nhat40_empcorr_ice(nhat20)[0]*n00,0,0,0,0], dtype=np.complex128)
    ODF_plots = (\
        {'nlm':nlm_constructor(-0.75), 'title':'', 'axloc':(0.205, 0.16), 'darr':arr(-90), 'lvlmax':0.4}, \
        {'nlm':nlm_constructor(+1.00), 'title':'', 'axloc':(0.577, 0.17), 'darr':arr(-90), 'lvlmax':0.4}, \
        {'nlm':nlm_constructor(+1.75), 'title':'', 'axloc':(0.737, 0.33), 'darr':arr(-90), 'lvlmax':0.4}, \
    )

    for ODF in ODF_plots:

        W = 0.125 # ax width
        axpos = [ODF['axloc'][0],ODF['axloc'][1], W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        nlm = ODF['nlm']
        lvlset = [np.linspace(0.0,ODF['lvlmax'],6), lambda x,p:'%.1f'%x]
        sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, showcb=False)

        # Arrow to ODF state
        n20_, n40_ = np.real(nlm[sf.I20])/norm, np.real(nlm[sf.I40])/norm
        ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                        xytext=(n20_+ODF['darr'][0]/norm, n40_+sc**2*ODF['darr'][1]/norm), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", linewidth=1.5, edgecolor='0.2', facecolor='0.2'),zorder=20)            

        axin.set_title(ODF['title'], fontsize=FS-3)

### Save figure

#plt.tight_layout()
plt.savefig('state-space-empcorr.png', transparent=1,  dpi=175)
plt.close()

