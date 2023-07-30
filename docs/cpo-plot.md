# Plotting

The orientation distribution function (ODF; normalized expansion series) can be plotted as follows:

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(4) 

import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cmasher as cmr
import cartopy.crs as ccrs

def plot_ODF(nlm, lm, ax, geo, cmap='Greys', cblabel='$n/N$ (ODF)', lvls=np.linspace(0.0,0.4,9), tickintvl=4, latres=60, plotAxes=False):
   
    # Discretize over S^2
    theta = np.linspace(0,   np.pi,   latres) # co-lat
    phi   = np.linspace(0, 2*np.pi, 2*latres) # lon
    phi, theta = np.meshgrid(phi, theta) # gridded 
    lon, colat = phi, theta
    lat = np.pi/2-colat
    _,nlm_len = lm.shape
    F = np.real(np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], phi,theta) for ii in np.arange(nlm_len) ], axis=0))
    F[F<0] = 0 # hide numerical/truncation errors
    
    # Plot    
    cmap = cmr.get_sub_cmap(cmap, 0.05, 1) # don't include pure white for visibility
    h = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='max', cmap=cmap, nchunk=5) # "nchunk" argument must be larger than 0 for constant-ODF (isotropy) to be plotted correctly.

    # Add grid lines
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

    # Colorbar
    cb = plt.colorbar(h, ax=ax, fraction=0.075, aspect=9,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])   
    cb.set_label(cblabel)
    cb.ax.xaxis.set_ticks(lvls, minor=True)
    
    if plotAxes:
        ax.plot([0],[90], marker=r'$z$', ms=9, c='tab:red',   transform=geo) # z axis
        ax.plot([90],[0], marker=r'$y$', ms=9, c='tab:blue',  transform=geo) # y axis
        ax.plot([0],[0],  marker=r'$x$', ms=9, c='tab:green', transform=geo) # x axis

    return h, cb   

### Plot ODF
    
# Make synthetic ODF
a2 = np.diag([0.0,0.5,0.5]) # an arbitrary a2
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # a2 contains information about the lowest-order harmonics l=2 only

# Setup figure
fig = plt.figure(figsize=(3,4))
inclination, rot = 45, +45 # view angle
prj, geo = ccrs.Orthographic(rot, 90-inclination), ccrs.Geodetic()
ax = plt.subplot(projection=prj)
ax.set_global() # show entire S^2

# Plot
h, cb = plot_ODF(nlm, lm, ax, geo, plotAxes=True)
fig.tight_layout()
plt.show()
```
