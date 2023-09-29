# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

import numpy as np
import scipy.special as sp
#import code # code.interact(local=locals())

from .discrete import cart2sph

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as mticker
import matplotlib.colors

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import cmasher as cmr
import cartopy.crs as ccrs

import warnings

### Discrete colors (colorbrewer)

c_red    = '#e31a1c'
c_lred   = '#fb9a99'
c_vlred  = '#fee0d2'
c_dred   = '#99000d'

c_blue   = '#1f78b4'
c_lblue  = '#a6cee3'
c_vlblue = '#deebf7'
c_dblue  = '#08519c'

c_green   = '#33a02c'
c_lgreen  = '#b2df8a'
c_vlgreen = '#e5f5e0'
c_dgreen  = '#006d2c'

c_lgray  = matplotlib.colors.to_hex('0.85')
c_gray   = matplotlib.colors.to_hex('gray')
c_dgray  = matplotlib.colors.to_hex('0.3')
c_vdgray = matplotlib.colors.to_hex('0.125')

c_orange = '#d94701'
c_yellow = '#ffff99'
c_purple = '#6a51a3'

c_brown  = '#b15928'
c_dbrown = '#8c510a'

### Common constants

fontsize_default = 12
isodist = 1/(4*np.pi) # value of ODF=n/N if isotropic and normalized

### Routines

def plotODF(nlm, lm, ax, \
        showcb=True, hidetruncerr=True, # options/flags \
        norm=True, latres=50, lonres=2*50, nchunk=5, # general \
        cmap='Greys', lvlset='iso-up', # colormap and lvls \ 
        cbfraction=0.075, cbaspect=9, cborientation='horizontal', cbpad=0.1, cblabel='$n/N$ (ODF)', cbtickintvl=4, # colorbar \
        kwargs_gridlines={'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}, # grid lines \
    ):
    
    """ 
    Plot distribution on S^2
    """

    ### Normalize 
    if norm is True:          N = np.sqrt(4*np.pi)*nlm[0] # default normalization; ODF=n/N
    elif type(norm) == float: N = norm # explicit normalization provided
    else:                     N = 1 # don't normalize
    nlm = np.divide(nlm, N)
    
    ### Discretize distribution on regular grid for contourf() 
    (F, lat, lon) = discretize(nlm, lm, latres, lonres)
    if hidetruncerr: F[F<0] = 0 # ignore numerical/truncation errors    

    ### Determine colormap
    if isinstance(cmap, str) and cmap == 'Greys' and \
       isinstance(lvlset, str) and lvlset == 'iso-up':
        cmap = cmr.get_sub_cmap('Greys', 0.125, 1) # don't include pure white to distinguish from white figure background
        cmap.set_under('w')
    
    ### Determine level set
    if isinstance(lvlset, str):
        if lvlset == 'iso-up':   
            lvls  = isodist*np.arange(1,9+0.1)
            lvlfmt = lambda x,pos:'$%i/(4\pi)$'%(x/isodist)
        elif lvlset == 'zero-up': 
            lvls  = isodist*np.arange(0,8+0.1)
            lvlfmt = lambda x,pos:'0' if x<1e-10 else '$%i/(4\pi)$'%(x/isodist)
        else: 
            raise ValueError('sfplt.plotODF(): Note sure what to do with passed lvlset; should be "iso-up", "zero-up" or list [lvls0, lvlmul, lvlfmt]')
    elif isinstance(lvlset, (list, tuple)) and len(lvlset) == 2:
        (lvls, lvlfmt) = lvlset # unpack and assume these are parse correctly
    else:
        print('lvlset=',lvlset)
        raise ValueError('sfplt.plotODF(): Note sure what to do with passed lvlset; should be "iso-up", "zero-up" or list [lvls0, lvlmul, lvlfmt]')

    ### Plot distribution
    # "nchunk" argument must be larger than 0 for isotropic distributions to be plotted correctly, else 0 is best choice.
    hodf = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, cmap=cmap, nchunk=nchunk, extend=('max' if lvls[0]<1e-10 else 'both')) 
    #ax.set_facecolor(color_bad) # Set different color for bad (masked) values (default white)

    ### Add grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

    ### Add colorbar
    if showcb:
        hcb = plt.colorbar(hodf, ax=ax, ticks=lvls[::cbtickintvl], fraction=cbfraction, aspect=cbaspect, orientation=cborientation, pad=cbpad)
        hcb.set_label(cblabel)
        cbax = hcb.ax.xaxis if cborientation == 'horizontal' else hcb.ax.yaxis
        cbax.set_ticks(lvls, minor=True)
        cbax.set_major_formatter(mticker.FuncFormatter(lvlfmt))
    else:
        hcb = None
        
    ### Return handles
    return hodf, hcb


def plotcoordaxes(ax, geo, axislabels='xi', color=c_dred, fontsize=None, negaxes=True):

    """
    Plot axis directions (labels)
    """

    if isinstance(axislabels, str):
        if axislabels   == 'xi':      lbls = ['$x$','$y$','$z$']
        elif axislabels == 'vuxi':    lbls = [r'$\vu{x}$',r'$\vu{y}$',r'$\vu{z}$']
        elif axislabels == 'vuei':    lbls = [r'$\vu{e}_1$',r'$\vu{e}_2$',r'$\vu{e}_3$']
        elif axislabels == 'vmi':     lbls = [r'$\vb{m}_1$',r'$\vb{m}_2$',r'$\vb{m}_3$']
        else: raise ValueError('sfplt.plotcoordinateaxes(): Note sure what to do with passed axislabels.')
    elif isinstance(axislabels, list) and len(axislabels) == 3:
        lbls = axislabels
    else:
        raise ValueError('sfplt.plotcoordinateaxes(): Note sure what to do with passed axislabels.')
    
    kwargs = dict(ha='center', va='center', transform=geo, color=color)
    if fontsize is not None: kwargs['fontsize'] = fontsize

    ax.text(0,  0, lbls[0], **kwargs)
    ax.text(90, 0, lbls[1], **kwargs)
    ax.text(0, 90, lbls[2], **kwargs)
            
    if negaxes:
        dphi = 180
        ax.text(dphi+0,  0, r'$-%s'%lbls[0][1:], **kwargs)
        ax.text(dphi+90, 0, r'$-%s'%lbls[1][1:], **kwargs)
            

def discretize(nlm, lm, latres, lonres):

    """
    Sample distribution on equispaced lat--lon grid
    """

    colat = np.linspace(0,   np.pi, latres) 
    lon   = np.linspace(0, 2*np.pi, lonres) 
    lon, colat = np.meshgrid(lon, colat) # gridded 
    lat = np.pi/2-colat

    nlmlen_from_nlm = len(nlm)
    nlmlen_from_lm = lm.shape[1]
    if nlmlen_from_nlm != nlmlen_from_lm: 
        nlmlen = np.amin([nlmlen_from_nlm,nlmlen_from_lm]) # pick smallest common range and continue (this is probably what the user wants)
        warnings.warn('sfplt.discretize(): dimensions of nlm (%i) and lm (%i) do not match, setting nlm_len=%i'%(nlmlen_from_nlm, nlmlen_from_lm, nlmlen))
    else:
        nlmlen = nlmlen_from_nlm
        
    F = np.real(np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], lon, colat) for ii in np.arange(nlmlen) ], axis=0))

    return (F, lat, lon)
    
    
def plotS2point(ax, v, *args, **kwargs):

    """
    Plot point on S2: wraps plt.plot()
    """

    lat, colat, lon = cart2sph(v, deg=True); 
    return ax.plot([lon],[lat], *args, **kwargs)


def getprojection(rotation=45, inclination=45):

    """
    Create orthographic projection class for above plotting
    """

    geo = ccrs.Geodetic()
    colatitude = 90-inclination
    prj = ccrs.Orthographic(rotation, colatitude)

    return geo, prj
    
    
def setfont_tex(fontfamily='serif', sansserif=['Times'], fontsize=fontsize_default):

    """
    Set latex-style fonts for math etc.
    """

    rc('font',**{'family':fontfamily, 'sans-serif':sansserif, 'size':fontsize})
    rc('text', usetex=True)
    rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage{txfonts} \usepackage{siunitx} \DeclareSIUnit\year{a}'

    return fontsize
    
 
def panellabel(ax, loc, txt, frameon=True, alpha=1.0, fontsize=fontsize_default, pad=0.3, ma='none', bbox=None, zorder=None):

    """
    Set panel label
    """

    at = AnchoredText(txt, loc=loc, prop=dict(size=fontsize), frameon=frameon, pad=pad, bbox_to_anchor=bbox, bbox_transform=ax.transAxes)
    at.patch.set_linewidth(0.7)
    ax.add_artist(at)


def plotFSE(ax, center, eigenvectors, eigenvalues, N=50, lw=0.5, ls='-', c='tab:green', scale=0.05):

    """
    Plot finite strain ellipse (FSE)
    """

    theta = np.linspace(0, 2*np.pi, N)
    eigenvectors = eigenvectors.T
    a = scale*np.sqrt(1/eigenvalues[0])
    b = scale*np.sqrt(1/eigenvalues[1])
    ellipse_points = a * np.cos(theta)[:, np.newaxis] * eigenvectors[:, 0] + b * np.sin(theta)[:, np.newaxis] * eigenvectors[:, 1]
    rotated_points = ellipse_points # Shift by center
    ax.plot(center[0] + rotated_points[:, 0], center[1] + rotated_points[:, 1],  lw=lw, ls=ls, c=c)
    
    
def plotparcel(ax, F, scale=1, axscale=1, elev=20, azim=35, \
                lw=1, facecolor='k', edgecolor='0.10',  \
                axislabels=True, colorax='k', fonttex=False, fontsize=fontsize_default, \
                drawinit=True, colorinit=c_dred, lwinit=1, \
                ):

    """
    Plot 3D parcel subject to deformation gradient F
    """

    ax.view_init(elev, 3*90-azim)
    
    ### Define vertices
    
    p_top   = np.array([ [0,0,1], [1,0,1], [1,1,1], [0,1,1] ]) * scale
    p_bot   = np.array([ [0,0,0], [1,0,0], [1,1,0], [0,1,0] ]) * scale
    p_left  = np.array([ [0,0,0], [0,0,1], [0,1,1], [0,1,0] ]) * scale
    p_right = np.array([ [1,0,0], [1,0,1], [1,1,1], [1,1,0] ]) * scale
    p_back  = np.array([ [0,1,0], [1,1,0], [1,1,1], [0,1,1] ]) * scale
    p_front = np.array([ [0,0,0], [1,0,0], [1,0,1], [0,0,1] ]) * scale

    ### Plot parcel sides
    
    args = (lw, '-', facecolor, edgecolor)
    at, ab, ah = 0.025, 0.10, 0.10
    _plotparcel_side(ax, p_top,   F, at, *args)
    _plotparcel_side(ax, p_bot,   F, ab, *args)
    _plotparcel_side(ax, p_front, F, ah, *args)
    _plotparcel_side(ax, p_left,  F, ah, *args)
    _plotparcel_side(ax, p_back,  F, ah, *args)
    _plotparcel_side(ax, p_right, F, ah, *args)

    ### Plot initial geometry

    if drawinit:
        kwargs = dict(ls='--', lw=lwinit, c=colorinit, zorder=10)
        points = ([0,0], [0,1*scale], [1*scale,1*scale], [1*scale,0])
        for p in points: ax.plot([p[0]]*2, [p[1]]*2, [0,scale*1], **kwargs)
        for p in points: ax.plot([0,scale*1], [p[0]]*2, [p[1]]*2, **kwargs)
        for p in points: ax.plot([p[0]]*2, [0,scale*1], [p[1]]*2, **kwargs)

    ### Axes 

    onespan = np.array([0,1])    
    ax.set_xlim(axscale*onespan) 
    ax.set_ylim(axscale*onespan)
    ax.set_zlim(axscale*onespan)
   
    ### x,y,z axes

    if axislabels:
        one = 0.5*scale*onespan
        kwargs = dict(color=colorax, fontsize=fontsize, va='center', ha='center', zorder=10)
        if fonttex: xilbl = [r"$\vu{x}$", r"$\vu{y}$", r"$\vu{z}$"]
        else:       xilbl = [r"$\bf{x}$", r"$\bf{y}$", r"$\bf{z}$"]
        ax.text(one.max() , 0, 0, xilbl[0], **kwargs)
        ax.text(0, one.max(), 0,  xilbl[1], **kwargs)
        ax.text(0, 0, one.max() , xilbl[2], **kwargs)

    ### Debug 
    
    if 0:
        args = (np.linspace(0,1,3),)
        ax.set_xticks(*args); 
        ax.set_yticks(*args); 
        ax.set_zticks(*args)
    else:
        ax.set_axis_off() # hide all parts of axes  
   

def _plotparcel_side(ax, pi, F, alpha, lw, ls, facecolor, edgecolor, zorder=5):

    qi = np.array([ np.matmul(F, pi[ii,]) for ii in range(4) ])
    x,y,z = qi[:,0],qi[:,1],qi[:,2]
    verts = np.array([list(zip(x,y,z))])
    coll = Poly3DCollection(verts)
    coll.set_edgecolor(edgecolor)
    coll.set_facecolor(facecolor)
    coll.set_linewidth(lw)
    coll.set_linestyle(ls)
    coll.set_alpha(alpha)
    coll.set_clip_on(False)
    coll.set_zorder(zorder)
    ax.add_collection3d(coll)
    
