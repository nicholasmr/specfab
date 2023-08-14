# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
import cmasher as cmr

from localheader import *

sys.path.append('../../')
import demolib as dl

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import constants as sfconst
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-0.5
FSANNO = FS-1.5

norm = 1/np.sqrt(4*np.pi)

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

### Setup figure

scale = 4.0
fig = plt.figure(figsize=(1.3*scale,0.9*scale))
plt.subplots_adjust(left=0.12, right=0.99, top=0.98, bottom=0.15)
ax = plt.gca()

xlims, ylims = [-1.45,2.65], [-1.4,3.5]
sc = np.diff(ylims)/np.diff(xlims)

### Determine valid subspace (valid eigenvalues)

print('Determining subspace of valid eigenvalues...')
isvalid, x,y = dl.nlm_isvalid_grid(xlims, ylims, RESX, RESY)

imdat = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
for ii in range(RESY):
    for jj in range(RESX):
        imdat[ii,jj,:-1] = matplotlib.colors.to_rgb('w') if isvalid[ii,jj] else matplotlib.colors.to_rgb(sfplt.c_lgray)

print('Plotting...')

# Plot valid/invalid subspaces and fabric-type shadings within valid subspace
imdat0 = np.ones((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)
kwargs = dict(aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)
im = ax.imshow(imdat0, **kwargs)
im = ax.imshow(imdat,  **kwargs)

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

### Plot ideal CPOs

ax.plot(n20_ideal, n40_ideal, '--', c='k', zorder=10, label=r'$\hat{n}_l^0 = Y_l^0(\theta,0)/\sqrt{4\pi}$')

### Misc

dl.plot_nlm_cases(ax, FSANNO, ms=8.5, dy0=0.07, show_circle=False)

plt.text(-1, 2.0, '{\\bf Unphysical}\n\\bf{eigenvalues}', color='0.3', ha='center', rotation=0, fontsize=FSANNO)

plt.sca(ax)
plt.xlabel(r'$\hat{n}_2^0$')
plt.ylabel(r'$\hat{n}_4^0$')

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'columnspacing': 0.5, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}
leg = plt.legend(loc=2, fontsize=FSLEG, frameon=False, ncol=2, **legkwargs); 

### Limits

plt.xlim(xlims)
plt.ylim(ylims)

### Plot ODF insets?

if 1:

    geo, prj = sfplt.getprojection(rotation=55+180, inclination=50)

    ### Modeled data
    
    arr = lambda ang: 0.125*np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
    n00 = norm
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

        W = 0.155 # ax width
        axpos = [ODF['axloc'][0],ODF['axloc'][1], W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        nlm = ODF['nlm']
        cmap = cmr.get_sub_cmap('Greys', 0.25, 1) # don't include pure white.
        cmap.set_under('w')
        lvlset = [np.linspace(ODF['lvlmin'], ODF['lvlmax'], 8), lambda x,p:'%.1f'%x]
        sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False)

        # Arrow to ODF state
        n20_, n40_ = np.real(nlm[3])/normfac, np.real(nlm[10])/normfac
        ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                        xytext=(n20_+ODF['darr'][0]/normfac, n40_+sc**2*ODF['darr'][1]/normfac), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", linewidth=1.5, edgecolor='0.2', facecolor='0.2'),zorder=20)            

        axin.set_title(ODF['title'], fontsize=FS)

### Save figure

plt.savefig('state-space-ideal.png', transparent=1,  dpi=175)

