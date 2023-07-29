# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

from localheader import *
sys.path.insert(0, '../../../demo') # for importing local specfabpy build (if available) and common python header
from header import * # contains matplotlib setup etc.
from specfabpy import specfabpy as sf

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

L=8
lm, nlm_len = sf.init(L) 

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
fig = plt.figure(figsize=(1.3*scale,0.9*scale))
plt.subplots_adjust(left=0.12, right=0.99, top=0.98, bottom=0.15)
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
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        if validregion[yii,xii]:
#            distnn = np.amin( np.sqrt(np.real((x_-x_corr[I])**2 + (1/sc*(y_-y_corr[I]))**2)) ) # distance from model-line points
#            var, expo = 1e-3, 6
#            imdat[yii,xii,-1] = 1 #np.exp(-distnn**expo/var) # set alpha depending on distance to model lines/points
            imdat[yii,xii,0:-1] = [1]*3 # OK
        else: 
            imdat[yii,xii,0:-1] = [0.85]*3 # bad 
            imdat[yii,xii,-1] = 1

print('done')

# Plot valid/invalid subspaces and fabric-type shadings within valid subspace
imdat0 = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
im = ax.imshow(imdat0, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
im = ax.imshow(imdat, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)

### Determine ideal boundary line

m = [0,0,1]
N = 100
thetavec = np.linspace(0, np.pi/2, N)
nlm_ideal = np.zeros((nlm_len, N))
n20_ideal, n40_ideal = np.zeros(N), np.zeros(N)
for ii, colat in enumerate(thetavec):
    nlm_ideal[:,ii] = np.real(sf.nlm_ideal(m, colat))
    n20_ideal[ii] = nlm_ideal[3,ii]/normfac
    n40_ideal[ii] = nlm_ideal[10,ii]/normfac

n20_unidir, n40_unidir = n20_ideal[0], n40_ideal[0]   # delta distributed
n20_planar, n40_planar = n20_ideal[-1], n40_ideal[-1] # x--y planar distributed

### End-member cases

# Isotropic state
ax.plot(0, 0, 'o', ms=mse, c='k', label=None, zorder=20)
dytext = 0.04/normfac
plt.text(0, 0+dytext, r'{\bf Isotropic}', color='k', ha='center', va='bottom', fontsize=FSANNO)

# Unidirectional
ax.plot(n20_unidir, n40_unidir, marker='o', ms=mse, ls='none', c=c_smax, label=None, zorder=20)
plt.text(n20_unidir-0.1, n40_unidir+dytext, '{\\bf Unidirectional}', color=c_smax, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# Planar 
ax.plot(n20_planar, n40_planar, marker='o', ms=mse, ls='none', c=c_girdle, label=None, zorder=20)
plt.text(n20_planar, n40_planar+dytext, '{\\bf Planar}', color=c_girdle, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# Shading labels
plt.text(-1, 2.0, '{\\bf Unphysical}\n\\bf{eigenvalues}', color='0.3', ha='center', rotation=0, fontsize=FSANNO)

### Plot ideal CPOs

ax.plot(n20_ideal, n40_ideal, '--', c='k', zorder=10, label=r'$\hat{n}_l^0 = Y_l^0(\theta,0)/\sqrt{4\pi}$')

### Aux

plt.sca(ax)
plt.xlabel(r'$\hat{n}_2^0$')
plt.ylabel(r'$\hat{n}_4^0$')

leg = plt.legend(loc=2, fontsize=FSLEG, frameon=False, ncol=2, **legkwargs); 

### Limits

plt.xlim(xlims)
plt.ylim(ylims)

### Plot ODF insets?

if 1:

    inclination = 50 # view angle
    rot0 = -90
    rot = -20 + rot0 
    prj = ccrs.Orthographic(rot, 90-inclination)
    geo = ccrs.Geodetic()     

    W = 0.155 # ax width
    tickintvl=1
    
    ### Modeled data
    
    arrmag = 0.125
    arr = lambda ang: arrmag*np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
    n00 = 1/np.sqrt(4*np.pi)
    I = [0, int(N/4*1.15), int(N/2), -1]
    lvlmax = 1.42
    ODF_plots = (\
        {'nlm':nlm_ideal[:,I[0]], 'title':r'$\theta=\SI{%i}{\degree}$'%(np.rad2deg(thetavec[I[0]])), 'axloc':(0.83, 0.56), 'darr':arr(-90), 'lvlmax':lvlmax}, \
        {'nlm':nlm_ideal[:,I[1]], 'title':r'$\theta=\SI{%i}{\degree}$'%(np.rad2deg(thetavec[I[1]])), 'axloc':(0.70, 0.17), 'darr':arr(-90), 'lvlmax':lvlmax}, \
        {'nlm':nlm_ideal[:,I[2]], 'title':r'$\theta=\SI{%i}{\degree}$'%(np.rad2deg(thetavec[I[2]])), 'axloc':(0.47, 0.31), 'darr':arr(+90), 'lvlmax':lvlmax}, \
        {'nlm':nlm_ideal[:,I[3]], 'title':r'$\theta=\SI{%i}{\degree}$'%(np.rad2deg(thetavec[I[3]])), 'axloc':(0.12, 0.24), 'darr':arr(-90), 'lvlmax':lvlmax}, \
    )

    for ODF in ODF_plots:

        axpos = [ODF['axloc'][0],ODF['axloc'][1], W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        nlm = ODF['nlm']
        lvls = np.linspace(0.0,ODF['lvlmax'],6)
    
        F, lon,lat = discretize_ODF(nlm, lm)
        F[F<0] = 0 # fix numerical/truncation errors
        h = axin.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend=('max' if lvls[0]==0.0 else 'both'), cmap='Greys', nchunk=5) # "nchunk" argument must be larger than 0 for constant-ODF (e.g. isotropy) is plotted correctly.

        # Arrow to ODF state
        n20_, n40_ = np.real(nlm[3])/normfac, np.real(nlm[10])/normfac
        ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                        xytext=(n20_+ODF['darr'][0]/normfac, n40_+sc**2*ODF['darr'][1]/normfac), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", linewidth=1.5, edgecolor='0.2', facecolor='0.2'),zorder=20)            

        # Add grid lines
        kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
        gl = axin.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
        gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

        axin.set_title(ODF['title'], fontsize=FS-1)

### Save figure

#plt.tight_layout()
plt.savefig('state-space-ideal.png', transparent=1,  dpi=175)
plt.close()

