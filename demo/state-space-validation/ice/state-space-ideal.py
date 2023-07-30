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

L = 8
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

print('Determining subspace of valid eigenvalues...')
xv, yv = np.meshgrid(x, y, indexing='xy')
validregion = np.reshape(sf.nlm_isvalid(xv.flatten(), yv.flatten()), (RESY, RESX))

### Determine subspace shadings

imdat = np.empty((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
imdat[:,:,0:-1] = [1.0]*3 # OK
            
print('Determining shading regions for different fabric types...')
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        if not validregion[yii,xii]: 
            imdat[yii,xii,0:-1] = [0.85]*3
            imdat[yii,xii,-1] = 1.0

print('Plotting...')

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
    nlm_ideal[:,ii] = np.real(sf.nlm_ideal(m, colat, L))
    n20_ideal[ii] = nlm_ideal[sf.I20,ii]/normfac
    n40_ideal[ii] = nlm_ideal[sf.I40,ii]/normfac

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
    lvlmax = 1.0
    lvlmin = 0.2
    title = lambda Ii: r'$\theta=\SI{%i}{\degree}$'%(np.rad2deg(thetavec[Ii]))
    ODF_plots = (\
        {'nlm':nlm_ideal[:,I[0]], 'title':title(I[0]), 'axloc':(0.83, 0.56), 'darr':arr(-90), 'lvlmax':lvlmax, 'lvlmin':lvlmin*1.8}, \
        {'nlm':nlm_ideal[:,I[1]], 'title':title(I[1]), 'axloc':(0.70, 0.17), 'darr':arr(-90), 'lvlmax':lvlmax, 'lvlmin':lvlmin}, \
        {'nlm':nlm_ideal[:,I[2]], 'title':title(I[2]), 'axloc':(0.47, 0.31), 'darr':arr(+90), 'lvlmax':lvlmax, 'lvlmin':lvlmin}, \
        {'nlm':nlm_ideal[:,I[3]], 'title':title(I[3]), 'axloc':(0.12, 0.24), 'darr':arr(-90), 'lvlmax':lvlmax, 'lvlmin':lvlmin}, \
    )

    for ODF in ODF_plots:

        axpos = [ODF['axloc'][0],ODF['axloc'][1], W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        nlm = ODF['nlm']
        lvls = np.linspace(ODF['lvlmin'], ODF['lvlmax'], 8)
        cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
        cmap.set_under('w')
    
        F, lon,lat = discretize_ODF(nlm, lm)
        F[F<0] = 0 # fix numerical/truncation errors
        h = axin.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend=('max' if lvls[0]==0.0 else 'both'), cmap=cmap, nchunk=5) 

        # Arrow to ODF state
        n20_, n40_ = np.real(nlm[3])/normfac, np.real(nlm[10])/normfac
        ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                        xytext=(n20_+ODF['darr'][0]/normfac, n40_+sc**2*ODF['darr'][1]/normfac), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", linewidth=1.5, edgecolor='0.2', facecolor='0.2'),zorder=20)            

        # Add grid lines
        kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
        gl = axin.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
        gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

        axin.set_title(ODF['title'], fontsize=FS)

### Save figure

#plt.tight_layout()
plt.savefig('state-space-ideal.png', transparent=1,  dpi=175)
plt.close()

