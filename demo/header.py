# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

"""
Shared routines used by several demo scripts
"""

import copy, sys, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc, colors
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker

import cmasher as cmr
import cartopy.crs as ccrs

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

### Setup matplotlib fonts etc.

FS = 12
FSSMALL = FS-1
rc('font',**{'family':'serif','sans-serif':['Times'],'size':FS})
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage{txfonts} \usepackage{siunitx} \DeclareSIUnit\year{a}'

### Default ODF plotting params

lvls_default = np.linspace(0.0,0.4,9)
tickintvl_default = 4

### Custom discrete colors (colorbrewer)

c_red    = '#e31a1c'
c_lred   = '#fb9a99'
c_dred   = '#99000d'
c_blue   = '#1f78b4'
c_lblue  = '#a6cee3'
c_green  = '#33a02c'
c_lgreen = '#b2df8a'
c_orange = '#d94701'
c_yellow = '#ffff99'
c_purple = '#6a51a3'
c_brown  = '#b15928'
c_gray   = 'gray'

"""
ODF plotting
"""

def plot_ODF(nlm, lm, ax=None, cmap='Greys', cblabel='$\psi$', rot0=-40, lvls=lvls_default, tickintvl=tickintvl_default, latres=60, nchunk=5):
    
    ### Plot ODF (harmonic distribution) on S^2 given vector of spectral coefficients (nlm)
    
    noaxis = ax is None # open new figure for plot?
    
    # Create new figure if needed
    if noaxis:
        size = 1.5
        plt.figure(figsize=(1.6*size,2*size))
        inclination = 45 # view angle
        rot = rot0 -90 # view angle
        prj = ccrs.Orthographic(rot, 90-inclination)
        geo = ccrs.Geodetic()
        ax = plt.subplot(projection=prj)
        ax.set_global()
    
    # Discretize distribution on regular grid so it can be plotted with contourf()
    F, lon,lat = discretize_ODF(nlm, lm, latres=latres)
    F[F<0] = 0 # fix numerical/truncation errors
    cmap = cmr.get_sub_cmap(cmap, 0.05, 1) # don't include pure white to distinguish from white figure background
    h = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend=('max' if lvls[0]==0.0 else 'both'), cmap=cmap, nchunk=nchunk) # "nchunk" argument must be larger than 0 for constant ODF (isotropy) to be plotted correctly.
    #ax.set_facecolor(color_bad) # Set different color for bad (masked) values (default white)

    # Add grid lines
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

    # Colorbar
    cb1 = plt.colorbar(h, ax=ax, fraction=0.075, aspect=9,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])   
    cb1.set_label(cblabel)
    cb1.ax.xaxis.set_ticks(lvls, minor=True)
    
    if noaxis: 
        plt.tight_layout()
        plt.show()
    
    return h


def discretize_ODF(nlm, lm, latres=60):

    ### Sample ODF (nlm) on regular grid needed for plotting

    theta = np.linspace(0,   np.pi,   latres) # CO-LAT 
    phi   = np.linspace(0, 2*np.pi, 2*latres) # LON
    phi, theta = np.meshgrid(phi, theta) # gridded 
    lon, colat = phi, theta # colat := theta 
    lat = np.pi/2-colat
    _,nlm_len = lm.shape
    F = np.real(np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], phi,theta) for ii in np.arange(nlm_len) ], axis=0))
    return (F, lon,lat)


"""
Modes of deformation: idealized strain histories with a constant velocity gradient

    F:      Deformation gradient tensor 
    ugrad:  Velocity gradient = d/dt(F) \dot F^-1
"""

class PureShear():

    # E.g. for ax='z': F = diag[ b^((1+r)/2), b^((1-r)/2), b^(-1) ]

    def __init__(self, t_c, r, ax='z'): 

        # e-folding time scale: 
        # if t_c > 0 ==> half-height time
        # if t_c < 0 ==> double-height time
        self.t_c = float(t_c)/np.log(2) 

        if not (0 <= r <= 1): 
            raise ValueError('Invalid bound, must be 0 <= r <= 1')

        # Exponents of diagonal entries of F
        rp, rm = (1+r)/2, (1-r)/2 # r=0 => uniaxial behaviour, else controls relative lengthening/shortening of other axes than "ax".
        if ax=='z': self.Fpow = [rp, rm, -1]
        if ax=='y': self.Fpow = [rp, -1, rm]
        if ax=='x': self.Fpow = [-1, rp, rm]

    # Deformation gradient
    def b(self, t): return np.exp(t/self.t_c) # scaling parameter (base of diagonal entries of F)
    def F(self, t): return np.diag(np.power(self.b(t),self.Fpow)) 
    
    # Strain tensor (cumulative)
    def strain(self, t): return 0.5*(self.F(t) + np.transpose(self.F(t))) - np.eye(3)

    # Strain-rate and spin tensors
    def D(self): return 1/self.t_c * np.diag(self.Fpow) # strain-rate
    def W(self): return np.diag([0,0,0])                # spin

    # Time taken to reach strain_zz strain for ax='z':
    #   strain_zz = F_zz - 1 = b^(-1) - 1 = exp(-t/t_e) - 1 => exp(-t/t_e) = strain_zz + 1
    def strainzz2time(self, strain_zz): return -self.t_c*np.log(strain_zz+1) 


class SimpleShear(): 

    # E.g. for plane='zx': F = [[1,0,tan(beta)],[0,1,0],[0,0,1]]

    def __init__(self, t_c, plane='zx'): 

        # Shear strain rate
        self.kappa = 1/t_c # t_c is the time taken to reach a shear strain of 1 (beta=45 deg.)

        # Shear plane
        if plane=='zx': self.Fshear = np.outer([1,0,0],[0,0,1])
        if plane=='yx': self.Fshear = np.outer([1,0,0],[0,1,0])

    # Deformation gradient
    def beta(self, time): return np.arctan(self.kappa*time) # beta = shear angle
    def F(self, t): return np.eye(3) + np.tan(self.beta(t))*self.Fshear
        
    # Strain tensor (cumulative)
    def strain(self, t): return 0.5*(self.F(t) + np.transpose(self.F(t))) - np.eye(3)
        
    # Strain-rate and spin tensors
    def ugrad(self): return self.kappa*self.Frel # velocity gradient
    def D(self): return 0.5*(self.ugrad() + self.ugrad().T) # strain-rate
    def W(self): return 0.5*(self.ugrad() - self.ugrad().T) # spin
    
    # Time taken to reach shear angle "beta"
    def beta2time(self, beta): return np.tan(beta)/self.kappa


"""
Plot 3D material parcel
"""

def plot_parcel(ax, xyz0, dzx,dzy,dyx, scale=1, elev=20, azim=70+180, color='k', colorax='k', lwax=1.25, lwside=1, FSAX=FSSMALL, plotaxlbls=False):
           
    # xyz0:         side lengths of parcel (x,y,z) 
    # dzx,dzy,dyx:  shear displacement of parcel surfaces
           
    x0,y0,z0 = xyz0 
    ax.view_init(elev, azim)
    kwargs = {'scale':scale, 'color':color, 'lw':lwside}
    
    ### Plot parcel sides
    
    x,y,z = [0,dyx,x0+dyx,x0], [0,y0,y0,0], [0,0,0,0] # bottom
    plot_side(ax, x,y,z, alpha=0.4, **kwargs)

    x,y,z = [0,dyx,dzx+dyx,dzx], [0,y0,y0+dzy,dzy], [0,0,z0,z0] # left
    plot_side(ax, x,y,z, alpha=0.3, **kwargs)
    
    x,y,z = [0,x0,x0+dzx,dzx,dzx], [0,0,dzy,dzy], [0,0,z0,z0] # back
    plot_side(ax, x,y,z, alpha=0.3, **kwargs)
    
    x,y,z = [dzx,dzx+dyx,x0+dzx+dyx,x0+dzx], [dzy,y0+dzy,y0+dzy,dzy], [z0,z0,z0,z0] # top
    plot_side(ax, x,y,z, alpha=0.1, **kwargs)
    
    x,y,z = [dyx,x0+dyx,x0+dzx+dyx,dzx+dyx], [y0,y0,y0+dzy,y0+dzy], [0,0,z0,z0] # front
    plot_side(ax, x,y,z, alpha=0.3, **kwargs)
    
    x,y,z = [x0,x0+dyx,x0+dzx+dyx,x0+dzx], [0,y0,y0+dzy,dzy], [0,0,z0,z0] # right
    plot_side(ax, x,y,z, alpha=0.3, **kwargs)
    
    ### x-y-z axes
    onespan = np.array([0,1])
    zero, one, doubleone = np.array([0,0]), np.array([0,1*scale]), np.array([1*scale,1*scale])
    xspan, yspan, zspan = one, one, one
    ax.plot(xspan,zero,zero, '-', lw=lwax, color=colorax, zorder=10)
    ax.plot(zero,yspan,zero, '-', lw=lwax, color=colorax, zorder=10)
    ax.plot(zero,zero,zspan, '-', lw=lwax, color=colorax, zorder=10)

    ### axes labels
    if plotaxlbls:
        axlenmul = 1.15
        ax.text(xspan.max()*axlenmul , 0, 0,                r"$\vu{x}$", color=colorax, fontsize=FSAX, zorder=10)
        ax.text(-xspan.max()*0.2, yspan.max()*axlenmul, 0,  r"$\vu{y}$", color=colorax, fontsize=FSAX, zorder=10)
        ax.text(0, 0, zspan.max()*axlenmul ,                r"$\vu{z}$", color=colorax, fontsize=FSAX, zorder=10)               
    
    ### For debugging 
#    ax.set_xticks([]) 
#    ax.set_yticks([]) 
#    ax.set_zticks([])
#    ax.set_xticks(np.linspace(0,1,5),minor=True)
#    ax.set_yticks(np.linspace(0,1,5),minor=True)
#    ax.set_zticks(np.linspace(0,1,5),minor=True)

    ax.set_xlim(onespan) 
    ax.set_ylim(onespan)
    ax.set_zlim(onespan)
    
    ax.set_axis_off() # hide all parts of axes  
   
    return


def plot_side(ax, x,y,z, alpha=0, scale=1, lw=1.0, color='k'):

    verts = scale*np.array([list(zip(x,y,z))])
    coll = Poly3DCollection(verts)
    coll.set_edgecolor('0.4')
    coll.set_facecolor(color)
    coll.set_linewidth(lw)
    coll.set_alpha(alpha)
    coll.set_clip_on(False)
    ax.add_collection3d(coll)
    
    
"""
Auxiliary routines
"""
    
def M_REG_custom(nu, expo, eps, sfobj):

    ### Custom regularization operator, M_REG

    nlm_len = int( (sfobj.Lcap+1)*(sfobj.Lcap+2)/2 )
    nlm_dummy = np.zeros((nlm_len), dtype=np.complex64)
    return -nu * np.linalg.norm(eps) * np.power(np.abs(sfobj.M_REG(nlm_dummy, eps)), expo)

    
def getPolarAngles(vec):
    x,y,z = vec
    phi   = np.rad2deg(np.arctan2(y,x))
    theta = 90 - np.rad2deg(np.arccos(z))
    return (theta, phi)
    

def writeSubplotLabel(ax, loc, txt, frameon=True, alpha=1.0, fontsize=FS, pad=0.3, ma='none', bbox=None):
    at = AnchoredText(txt, loc=loc, prop=dict(size=fontsize), frameon=frameon, pad=pad, bbox_to_anchor=bbox, bbox_transform=ax.transAxes)
    at.patch.set_linewidth(0.7)
    ax.add_artist(at)

