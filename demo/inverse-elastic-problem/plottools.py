# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Library providing some basic plotting tools used in several scripts
"""

import sys
sys.path.insert(0, '..')
from header import * # parent plotting libraries 
from inverseproblem import cart2sph

import matplotlib.pyplot as plt

def plot_ODF(nlm, lm, ax=None, cmap='Greys', cblabel='$\psi$', rot0=-40, \
                lvls=lvls_default, tickintvl=tickintvl_default, \
                cbaspect=7, cbfrac=0.1, cborientation='horizontal'):

    # Plot an ODF prescribed in terms of "nlm" (similar to the routine in parent library, but this overrides it)
        
    F, lon,lat = discretize_ODF(nlm, lm)
    F[F<0] = 0 # fix numerical/truncation errors
    cmap = cmr.get_sub_cmap(cmap, 0.05, 1) # don't include pure white.
    h = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend=('max' if lvls[0]==0.0 else 'both'), cmap=cmap, nchunk=5) # "nchunk" argument must be larger than 0 for isotropic ODFs to be plotted correctly

    # Add grid lines
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

    # Colorbar
    cb1 = plt.colorbar(h, ax=ax, fraction=cbfrac, aspect=cbaspect,  orientation=cborientation, pad=0.1, ticks=lvls[::tickintvl])   
    cb1.set_label(cblabel)
    if cborientation == 'horizontal': cb1.ax.xaxis.set_ticks(lvls, minor=True)
    else:                             cb1.ax.yaxis.set_ticks(lvls, minor=True)
    
    return h

def get_deg(lat_or_colat, lon):
    return (np.rad2deg(lat_or_colat), np.rad2deg(lon))

def plot_unitaxes(ax, geo, colorxi='#662506'):
    # Plot x,y,z axes on S^2 
    kwargs = {'deg':True, 'colat':True}
    vx,vy,vz = cart2sph([1,0,0], **kwargs), cart2sph([0,1,0], **kwargs), cart2sph([0,0,1], **kwargs)
    vym = cart2sph([0,-1,0], **kwargs)
    ax.text(vx[1], vx[0], r'$\vu{x}$', color=colorxi, horizontalalignment='center', transform=geo)
    ax.text(vy[1], vy[0], r'$\vu{y}$', color=colorxi, horizontalalignment='center', transform=geo)
    ax.text(vym[1], vym[0], r'$-\vu{y}$', color=colorxi, horizontalalignment='center', transform=geo)
    ax.text(vz[1], vz[0], r'$\vu{z}$', color=colorxi, horizontalalignment='center', transform=geo)
    
