#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021

import copy, sys, os, code # code.interact(local=locals())
import numpy as np

from scipy import interpolate
import scipy.special as sp

from specfabpy import specfabpy as sf # requires the spectral fabric module to be compiled!

s2yr   = 3.16887646e-8
yr2s   = 31556926    
yr2kyr = 1e-3 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import rcParams, rc, patches
import matplotlib.ticker as mticker
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs

FS = 12
rc('font',**{'family':'serif','sans-serif':['Times'],'size':FS})
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage{txfonts} \usepackage{siunitx}'

# plot colors (colorbrewer)
cb = '#1f78b4'
cr = '#e31a1c'
cg = '#33a02c'

colax = '#99000d'

class PureShear():

    def __init__(self, t_c, r, ax='z'): 

        # e-folding time scale: 
        # if t_c > 0 ==> half-height time
        # if t_c < 0 ==> double-height time
        self.t_c = float(t_c)/np.log(2) 

        if ax=='z': self.Fpow = [(1+r)/2., (1-r)/2., -1]
        if ax=='y': self.Fpow = [(1+r)/2., -1, (1-r)/2.]
        if ax=='x': self.Fpow = [-1, (1+r)/2., (1-r)/2.]

    def lam(self, t):    return np.exp(t/self.t_c) # lambda(t)
    def F(self, t):      return np.diag(np.power(self.lam(t),self.Fpow))
    
    def strain(self, t): return 0.5*( self.F(t) + np.transpose(self.F(t)) ) - np.diag([1,1,1])
    #def epszz2time(self,epszz): return -self.t_e*np.log(epszz+1) # time it takes to reach "eps_zz" strain with t_e char. timescale

    # Note that F is constructed such that W and eps are time-independant.
    def D(self): return 1/self.t_c * np.diag(self.Fpow)
    def W(self): return np.diag([0,0,0])
    
def get_mycmap(vmax, Ncolors):

    vmin = 0
    kwargs = {"vmin":vmin, "vmax": vmax}
    lvls = np.linspace(vmin,vmax,Ncolors+1); # print lvls

    import brewer2mpl
    from matplotlib import rcParams
    import matplotlib as mpl 

    dN = 1
    RAWMAP = brewer2mpl.get_map('Greys', 'Sequential', Ncolors+dN, reverse=False).hex_colors
    cmap = mpl.colors.ListedColormap(RAWMAP[:-dN])
    cmap.set_over(RAWMAP[-dN])
    
    return (lvls,cmap,kwargs)

def plot_ODF(ax, nlm, lm, geo=None, rot0=0, cmap='Greys', tickintvl=2, cblbl=r'ODF', title=''):
    
    F, lon,lat = discretize_fabric(nlm, lm)
    F[np.where(F<0)] = 0 # trancation error might give rise to small negative values in the ODF, don't show these neglieble errors.
    if 1:
        k = 13; q = 9
        F[q:(q+1),k:(k+1)] = np.amin(F[1:2,k:(k+1)])*10

#    lvls = np.linspace(0.1,0.8,8) # Contour lvls
    (lvls,cmap,kwargs) = get_mycmap(0.2,4)
    
    hdistr = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='max', cmap=cmap, **kwargs)
    kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))
    cb1 = plt.colorbar(hdistr, ax=ax, fraction=0.065, aspect=10,  orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])   
    cb1.set_label(cblbl)
    cb1.ax.xaxis.set_ticks(lvls, minor=True)
    ax.set_title(title, fontsize=FS)
    
    if geo is not None:
        ax.plot([0],[90],'.', c=colax, transform=geo)
        ax.plot([rot0-90],[0],'.', c=colax, transform=geo)
        ax.plot([rot0],[0],'.', c=colax, transform=geo)
        ax.text(rot0-80, 70, r'$\vu{z}$', c=colax, horizontalalignment='left', transform=geo)
        ax.text(rot0-90+10, -5, r'$\vu{x}$', c=colax, horizontalalignment='left', transform=geo)
        ax.text(rot0-15, -5, r'$\vu{y}$', c=colax, horizontalalignment='left', transform=geo)
    
    return hdistr

def discretize_fabric(nlm, lm, latres=50):

    theta = np.linspace(0,   np.pi,   latres) # CO-LAT 
    phi   = np.linspace(0, 2*np.pi, 2*latres) # LON
    phi, theta = np.meshgrid(phi, theta) # gridded 
    lon, colat = phi, theta
    lat = np.pi/2-colat
    _,nlm_len = lm.shape
    F = np.real(np.sum([ nlm[ii]*sp.sph_harm(lm[1][ii], lm[0][ii], phi,theta) for ii in np.arange(nlm_len) ], axis=0))
    
    return (F, lon,lat)

#-------------------
# SynthFab
#-------------------

class SyntheticFabric():
    
    def __init__(self, L=16, nu=3e-3): 

        self.L  = L # spectral truncation 
        self.nu0 = nu # Regularization strength (calibrated for L and strain-rate)
               
        self.nlm_len = sf.init(self.L) 
        self.lm = sf.get_lm(self.nlm_len)
        
    def intp_nlm(self, nlm_list, depth):
        
        nlm_list_intp = np.zeros((self.N,self.nlm_len), dtype=np.complex128) 
        z_intp = np.linspace(0, depth[-1], self.N)
        
        for lmii in np.arange(self.nlm_len):
        
            intp = interpolate.interp1d(depth, nlm_list[:,lmii])
            nlm_list_intp[:,lmii] = intp(z_intp)
        
        return (nlm_list_intp, -z_intp)
    
    def plot_side(self, ax, x,y,z, alpha=0, scale=1, lw=1.25, color='k'):
        verts = scale*np.array([list(zip(x,y,z))])
        coll = Poly3DCollection(verts)
        coll.set_edgecolor('0.4')
        coll.set_facecolor(color)
        coll.set_linewidth(lw)
        coll.set_alpha(alpha)
        coll.set_clip_on(False)
        ax.add_collection3d(coll)
    
    def plot_parcel(self, ax, xyz0, dzx,dzy,dyx, scale=0.7, color='k', plotaxlbls=False):
               
        x0,y0,z0   = xyz0 
        ax.view_init(20, +70+180)
        
        x,y,z = [0,dyx,x0+dyx,x0], [0,y0,y0,0], [0,0,0,0] # bottom
        self.plot_side(ax, x,y,z, alpha=0.4, scale=scale, color=color)

        x,y,z = [0,dyx,dzx+dyx,dzx], [0,y0,y0+dzy,dzy], [0,0,z0,z0] # left
        self.plot_side(ax, x,y,z, alpha=0.3, scale=scale, color=color)
        
        x,y,z = [0,x0,x0+dzx,dzx,dzx], [0,0,dzy,dzy], [0,0,z0,z0] # back
        self.plot_side(ax, x,y,z, alpha=0.3, scale=scale, color=color)
        
        x,y,z = [dzx,dzx+dyx,x0+dzx+dyx,x0+dzx], [dzy,y0+dzy,y0+dzy,dzy], [z0,z0,z0,z0] # top
        self.plot_side(ax, x,y,z, alpha=0.1, scale=scale, color=color)
        
        x,y,z = [dyx,x0+dyx,x0+dzx+dyx,dzx+dyx], [y0,y0,y0+dzy,y0+dzy], [0,0,z0,z0] # front
        self.plot_side(ax, x,y,z, alpha=0.3, scale=scale, color=color)
        
        x,y,z = [x0,x0+dyx,x0+dzx+dyx,x0+dzx], [0,y0,y0+dzy,dzy], [0,0,z0,z0] # right
        self.plot_side(ax, x,y,z, alpha=0.3, scale=scale, color=color)
        
        #ax.scatter([0],[0],[0], 'o', s=[50], color='k')
        lw=1.25
        onespan = np.array([0,1])
        zero, one = np.array([0,0]), np.array([0,1*scale])
        xspan, yspan, zspan = one, one, one
        ax.plot(xspan,zero,zero, '-', lw=lw, color=colax, zorder=10)
        ax.plot(zero,yspan,zero, '-', lw=lw, color=colax, zorder=10)
        ax.plot(zero,zero,zspan, '-', lw=lw, color=colax, zorder=10)
        if plotaxlbls:
            axlenmul = 1.15
            ax.text(xspan.max()*axlenmul , 0, 0, r"$\vu{x}$", color=colax, zorder=10)
            ax.text(-xspan.max()*0.2, yspan.max()*axlenmul, 0, r"$\vu{y}$", color=colax, zorder=10)
            ax.text(0, 0, zspan.max()*axlenmul , r"$\vu{z}$", color=colax, zorder=10)               
        
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); 
        ax.set_xticks(np.linspace(0,1,5),minor=True)
        ax.set_yticks(np.linspace(0,1,5),minor=True)
        ax.set_zticks(np.linspace(0,1,5),minor=True)
        ax.set_xlim(onespan); ax.set_ylim(onespan); ax.set_zlim(onespan)
        
        ax.set_axis_off()
        
       
        return
    
    ####################################

    def makeProfile(self, dt=25):
    
        deformExpr1 = {'ax':'z', 't_c':+1000, 'r':0, 'Gamma0':0}
        deformExpr = deformExpr1
        
        t_end = [3000*yr2s, 1500*yr2s]
        dt    *= yr2s 
        
        Nt0 = int(t_end[0]/dt)# Total number of integration steps taken
        Nt1 = int(t_end[1]/dt)# Total number of integration steps taken
        Nt = Nt0 + Nt1
        
        ### Construct strain-rate and spin tensor histories
        
        D = np.zeros((Nt+1,3,3)) # strain-rate
        W = np.zeros((Nt+1,3,3)) # spin
        
        # Parcel deformation for plotting
        xyz0_init = [1,1,1]
#        xyz0 = np.repeat(xyz0_init, Nt)
        
        xyz0     = np.ones((Nt+1,3))
        xyz0_alt = np.ones((Nt+1,3))
        
        strain_zz = np.zeros(Nt+1)
        
        # Determine D and W for fabric evolution
        for ii in np.arange(Nt0+1):
            t = ii*dt
            PS = PureShear(+deformExpr['t_c']*yr2s, deformExpr['r'], ax=deformExpr['ax'])
            strain_zz[Nt0-ii] = PS.strain(t)[-1,-1]
            xyz0[Nt0-ii,:] = np.matmul(PS.F(t), xyz0_init)
            W[Nt0-ii,:,:], D[Nt0-ii,:,:] = PS.W(), PS.D()            

        for ii in np.arange(Nt1+1):
            t = ii*dt
            PS = PureShear(-deformExpr['t_c']*yr2s, deformExpr['r'], ax=deformExpr['ax'])
            strain_zz[Nt0+ii] = PS.strain(t)[-1,-1]
            xyz0[Nt0+ii,:] = np.matmul(PS.F(t), xyz0_init)
            W[Nt0+ii,:,:], D[Nt0+ii,:,:] = PS.W(), PS.D()            

        print(Nt, Nt1,Nt0, strain_zz)

        ### Fabric evolution 
        
        nlm_list      = np.zeros((Nt+1,self.nlm_len), dtype=np.complex128) # The expansion coefficients
        nlm_list[:,0] = 1/np.sqrt(4*np.pi) # Normalized such that N(t=0) = 1
        eigvals = np.zeros((Nt+1,3))
        
        Eij = np.ones((Nt+1,3,3))
        e1,e2,e3 = [1, 0, 0], [0, 1, 0], [0, 0, 1]
        Ecc,Eca,alpha,nprime = sf.ecc_opt_lin, sf.eca_opt_lin, sf.alpha_opt_lin, 1
        #Ecc,Eca,alpha,nprime = 1, 1e4, 0, 1
        
        nlm_prev = nlm_list[Nt,:].copy()
        for ii in np.arange(Nt0+1):
            tt = Nt0-ii
            D_, W_ = D[tt,:,:], W[tt,:,:]
            dndt = sf.nu(self.nu0, D_) * sf.dndt_reg(nlm_prev) # Regularization
            dndt += sf.dndt_latrot(nlm_prev, D_,W_) # Lattice rotation
            nlm_list[tt,:] = nlm_prev + dt * np.matmul(dndt, nlm_prev)
            nlm_prev = nlm_list[tt,:]
            Eij[tt,:,:] = sf.enhfac_eiej(nlm_prev, e1,e2,e3, Ecc,Eca,alpha,nprime) 
            _,_,_, eigvals[tt,:] = sf.frame(nlm_prev, 'e')
#            print(eigvals)
        
        nlm_prev = nlm_list[Nt,:].copy()
        for ii in np.arange(Nt1+1):
            tt = Nt0+ii
            D_, W_ = D[tt,:,:], W[tt,:,:]
            dndt = sf.nu(self.nu0, D_) * sf.dndt_reg(nlm_prev) # Regularization
            dndt += sf.dndt_latrot(nlm_prev, D_,W_) # Lattice rotation
            nlm_list[tt,:] = nlm_prev + dt * np.matmul(dndt, nlm_prev)
            nlm_prev = nlm_list[tt,:]
            Eij[tt,:,:] = sf.enhfac_eiej(nlm_prev, e1,e2,e3, Ecc,Eca,alpha,nprime) 
            _,_,_, eigvals[tt,:] = sf.frame(nlm_prev, 'e')
        
        ### Plot
        
        steps = [0, int(Nt/2),Nt]
        #steps = np.arange(0,Nt,10)
        steps = [0,Nt0,Nt]
        m = 1
        steps = np.hstack((np.arange(0,Nt0,m),np.arange(Nt0,Nt+1,m))) 
        
        for ii in steps:
            
            #scale = 0.4
            scale = 0.4
            fig = plt.figure(figsize=(14*scale,10*scale))
            
            gs_master = gridspec.GridSpec(2, 2, width_ratios=[0.8,1], height_ratios=[1,0.8])
            gs_master.update(top=0.91, bottom=0.13, left=0.12, right=1-0.15, hspace=-0.05, wspace=0.1)
            gs = gs_master
 
            axE = fig.add_subplot(gs_master[:, 0])
            #
            axParcel = fig.add_subplot(gs[0, 1], projection='3d')
            #
            inclination = 50 # view angle
            rot0 = 1 * -90
            rot = 1*-30 + rot0 #*1.4 # view angle
            prj = ccrs.Orthographic(rot, 90-inclination)
            geo = ccrs.Geodetic()            
            axODF = fig.add_subplot(gs[1, 1], projection=prj)
            axODF.set_global()
            
            ##
 
            lw = 1.5
            axE.plot(Eij[:,2,2],strain_zz, '--k',  lw=lw, label=r'$E_{zz}$')
            axE.plot(Eij[:,0,2],strain_zz, '-k', lw=lw, label=r'$E_{xz}$')

            ylims = [-1,2]
            dy=0.5
            axE.set_yticks(np.arange(ylims[0],ylims[1]+dy, dy))
            axE.set_yticks(np.arange(ylims[0],ylims[1]+dy, dy/2), minor=True)
            axE.set_ylim(ylims)

            xlims = [0,2.5]
            axE.plot(xlims,strain_zz[[ii,ii]], ':', c=cg, lw=lw)
            dx=1
            axE.set_xticks(np.arange(xlims[0],xlims[1]+dx,dx))
            axE.set_xticks(np.arange(xlims[0],xlims[1]+dx,dx/4), minor=True)
            axE.set_xlim(xlims)
            
            axE.set_ylabel('$\epsilon_{zz}$')
            axE.set_xlabel('$E_{ij}$')
            axE.set_title('Directional enhancement factors', pad=12, fontsize=FS)
            #
            legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':0.9, 'ncol':1, 'handlelength':1.34, 'labelspacing':0.3}
            hleg = axE.legend(loc=1, **legkwargs)
            hleg.get_frame().set_linewidth(0.7);
            
            ###

            axE.text(1.95,1.05, r'Rathmann et al. (2021)', transform=axE.transAxes, fontsize=FS-1, horizontalalignment='left')
            
            y0,dy = 0.20, 0.06
            axE.text(2.2,y0+2*dy, r'$\lambda_1 = %.2f$'%(eigvals[ii,0]), transform=axE.transAxes, fontsize=FS-1, ha='left')
            axE.text(2.2,y0+1*dy, r'$\lambda_2 = %.2f$'%(eigvals[ii,1]), transform=axE.transAxes, fontsize=FS-1, ha='left')
            axE.text(2.2,y0+0*dy, r'$\lambda_3 = %.2f$'%(eigvals[ii,2]), transform=axE.transAxes, fontsize=FS-1, ha='left')
            
            self.plot_parcel(axParcel, xyz0[ii,:], 0,0,0, color='0.5', plotaxlbls=True)
            plot_ODF(axODF, nlm_list[ii,:], self.lm, cmap='Greys', geo=geo, rot0=rot0)        
        
            ###

            print('Saving %i'%(ii))        
            if ii <= Nt0:
                fout1 = 'cube_crush_frames/frame%04i.png'%(Nt0-ii)
                fout2 = 'cube_crush_frames/frame%04i.png'%(Nt0+ii)
                plt.savefig(fout1, dpi=300)
                os.system('cp %s %s'%(fout1,fout2))

            if ii > Nt0:
                fout1 = 'cube_crush_frames/frame%04i.png'%(Nt0 +ii)
                fout2 = 'cube_crush_frames/frame%04i.png'%(2*Nt - (ii-Nt0))
                plt.savefig(fout1, dpi=300)
                print('2*Nt=%i :: %i :: %i'%(2*Nt, Nt0 +ii, 2*Nt - (ii-Nt0)))
                os.system('cp %s %s'%(fout1,fout2))

            plt.close()
    
###

os.system('mkdir -p cube_crush_frames')    
synfab = SyntheticFabric()
synfab.makeProfile()    
os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 3 -i cube_crush_frames/frame%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p cube_crush_ani.avi')
    
       
