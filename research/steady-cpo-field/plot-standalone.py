#!/usr/bin/python3
# Nicholas Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Standalone plotter for results of steady SSA fabric solver
"""

import numpy as np
import sys, code, pickle
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

ms2myr = 3.17098e+8
m2km = 1e-3
mapscale = m2km # x and y axis scale

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True) 

"""
Variable definitions
--------------------
    ux[:], uy[:]   : x,y horizontal velocity components (m/s)
    umag[:]        : velocity vector magnitude (m/s)
    epsE[:]        : effective strain rate (1/s)
    S[:],B[:],H[:] : surface, bottom and height of ice column
    mask[:]        : domain mask, values 1-5: ocean (1), ice_free_land (2) grounded_ice (3) floating_ice (4)
    mi[:,i]        : i-th a^(2) eigenvector, sorted so that i=0,1 correspond to the largest and second-largest eigenvalue **in the model plane**, respectively, and i=2 is the out-of-plane eigenvector
    lami[:,i]      : i-th a^(2) eigenvalue, sorted similarly to mi
    E[:]           : equivelent isotropic enhancement factor, estimated using CAFFE
"""

### Load solution

if len(sys.argv) != 3:
    print('usage: %s inputs.pkl solution.pkl'%(sys.argv[0]))
    sys.exit(1)
    
finp = sys.argv[1]
fsol = sys.argv[2]
coords,cells, ux,uy,umag,epsE, S,B,H,mask = pickle.load(open(finp, 'rb'))
coords,cells, mi,lami,E_CAFFE = pickle.load(open(fsol, 'rb'))

coords = (mapscale*coords[0], mapscale*coords[1]) # rescale coordinates
triang = tri.Triangulation(*coords, triangles=cells) # triangulation object for FEM mesh

### Plot

xlims = [np.amin(coords[0]),np.amax(coords[0])]
ylims = [np.amin(coords[1]),np.amax(coords[1])]

def newfig(floating=True, mesh=False, figsize=(5,5)):
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    legh, legt = [], []
    if mesh: 
        ax.triplot(triang, lw=0.075, color='0.5', alpha=0.8, zorder=10)
    if floating: 
        ax.tricontour(triang, mask==3, [0.5, 1.5], colors=['limegreen',], linewidths=2, zorder=11)
        legh.append(Line2D([0], [0], color='limegreen', lw=2))
        legt.append('Floating')
    ax.legend(legh, legt, ncol=2, loc=1, bbox_to_anchor=(1,1.15), handlelength=1.2, fancybox=False, frameon=False)
    ax.axis('square')
    ax.set_xlabel(r'$x$ (km)')
    ax.set_ylabel(r'$y$ (km)')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.add_patch(plt.Rectangle((xlims[0],ylims[0]), np.diff(xlims)[0], np.diff(ylims)[0], color='0.85'))
    return (fig, ax)
    
def newcax(ax): return make_axes_locatable(ax).append_axes("right", size="4%", pad=0.13)

kw_save = dict(dpi=150, pad_inches=0.1, bbox_inches='tight')

### Velocities

fig, ax = newfig()
lvls = np.logspace(0.5, 3.5, 13)
cs = ax.tricontourf(triang, ms2myr*umag, levels=lvls, norm=colors.LogNorm(vmin=lvls[0], vmax=lvls[-1]), extend='both', cmap='inferno')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label('$u$ (m/yr)')
plt.savefig('umag.png', **kw_save)

### Effective strain rate

fig, ax = newfig()
lvls = np.arange(0, 50+.01, 5)
cs = ax.tricontourf(triang, 1e3*ms2myr*epsE, levels=lvls, extend='max', cmap='viridis')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$\dot{\epsilon}_{e}$ (1/yr)')
plt.savefig('epsE.png', **kw_save)

### Delta lambda (horizontal eigenvalue difference)

fig, ax = newfig()
lvls = np.arange(0, 0.8+.01, 0.1)
dlam = abs(lami[:,0] - lami[:,1]) # eigenvalues 1 and 2 are the largest and smallest in-model-plane eigenvalues
cs = ax.tricontourf(triang, dlam, levels=lvls, extend='max', cmap='Blues')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$\Delta\lambda$')

# Quiver principal horizontal eigenvector
meshpts = (coords[0], coords[1])
xv, yv = np.linspace(*xlims, 15)[1:-1], np.linspace(*ylims, 15)[1:-1]
x, y = np.meshgrid(xv, yv, indexing='xy')
m1 = mi[:,0,:] # principal horizontal eigenvector
m1x = griddata(meshpts, m1[:,0].flatten(), (x, y), method='linear', fill_value=np.nan)
m1y = griddata(meshpts, m1[:,2].flatten(), (x, y), method='linear', fill_value=np.nan) # y coordinate is index 2 (z coordinate) since problem is in xz plane
renorm = np.sqrt(m1x**2+m1y**2)
m1x, m1y = np.divide(m1x, renorm), np.divide(m1y, renorm)
hq = ax.quiver(x, y, +m1x, +m1y, color='tab:red', scale=60)
hq = ax.quiver(x, y, -m1x, -m1y, color='tab:red', scale=60)
ax.quiverkey(hq, 0.05, 1.05, 3, r'${\bf m}_1$', labelpos='E')

plt.savefig('dlam.png', **kw_save)

### lambda_z (vertical eigenvalue)

fig, ax = newfig()
lvls = np.arange(0, 0.8+.01, 0.1)
lamz = lami[:,2] # eigenvalue 3 is the out-of-model-plane (z) eigenvalue
cs = ax.tricontourf(triang, lamz, levels=lvls, extend='max', cmap='RdPu')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$\lambda_z$')
plt.savefig('lamz.png', **kw_save)

### E (CAFFE)

fig, ax = newfig()
lvls = np.logspace(-1, 1, 17)
cs = ax.tricontourf(triang, E_CAFFE, levels=lvls, norm=colors.LogNorm(vmin=lvls[0], vmax=lvls[-1]), extend='both', cmap='PuOr_r')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$E$', labelpad=-2)
plt.savefig('E.png', **kw_save)

