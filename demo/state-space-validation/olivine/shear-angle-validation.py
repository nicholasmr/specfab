# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2024

import copy, sys, code # code.interact(local=locals())

import numpy as np
import pickle

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmasher as cmr
import cartopy.crs as ccrs

sys.path.insert(0, '..')
from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)


from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import integrator as sfint
from specfabpy import common as sfcom
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-1

#--------------------
# Run options
#--------------------

Ll = 6  # L low
Lh = 14 # L high

Nt = 250

#--------------------
# Model integratation
#--------------------

### specfab runs

mod, target = dict(type='ss', plane=1, T=1), np.deg2rad(84)

kwargs   = dict(Nt=Nt, iota=-1, zeta=0, nu=1, Gamma0=None, Lambda=None)
kwargs_D = dict(Nt=Nt, iota=-1, zeta=0, nu=1, Gamma0=None, Lambda=2.5e-2)

lm_Ll, nlm_len_Ll = sf.init(Ll) 
sf_nlm_Ll,  sf_F_Ll,  sf_time_Ll,  sf_ugrad_Ll  = sfint.lagrangianparcel(sf, mod, target, **kwargs) 
sfd_nlm_Ll, sfd_F_Ll, sfd_time_Ll, sfd_ugrad_Ll = sfint.lagrangianparcel(sf, mod, target, **kwargs_D)

lm_Lh, nlm_len_Lh = sf.init(Lh) 
sf_nlm_Lh,  sf_F_Lh,  sf_time_Lh,  sf_ugrad_Lh  = sfint.lagrangianparcel(sf, mod, target, **kwargs) 
sfd_nlm_Lh, sfd_F_Lh, sfd_time_Lh, sfd_ugrad_Lh = sfint.lagrangianparcel(sf, mod, target, **kwargs_D) 

### D-Rex runs

drex_nlm, drex_F = load_drex_run('simpleShear', 1)

### Skemer et al. data

ZhangKarato1995_1200C_A_strain = [17,30,45,65,110,]
ZhangKarato1995_1200C_A_angle  = [43,37,38,24,20,]

ZhangKarato1995_1300C_A_strain = [11,7,65,58,100,115,150,]
ZhangKarato1995_1300C_A_angle  = [36,28,18,10,5,10,0,]

Warren2008_A_strain = [0,65,118,131,258,386,386,525,168,]
Warren2008_A_angle  = [62,37,49,61,4,11,0,1,33,] 

Skemer2011_A_strain = [120,180,350,]
Skemer2011_A_angle  = [55,53,45,]

Webber2010_A_strain = [0,25,100,130,168,330,330,]
Webber2010_A_angle  = [55,35,47,29,37,47,45,]

Jung2006_A_strain = [120,]
Jung2006_A_angle  = [26,]

#Jung2006_B_strain = [400,]
#Jung2006_B_angle  = [15,]

#Jung2006_C_strain = [100,110,120,400,]
#Jung2006_C_angle  = [27,25,16,36,]

Hansen2014_A_strain = [500,590,680,680,760,820,870,880,1020,1060,1090,1420,1870,]
Hansen2014_A_angle  = [5.4,0.6,-5.5,-4.30000000000001,-8.40000000000001,-1.90000000000001,-0.900000000000006,-8.90000000000001,2.8,-2,-1.5,-4.09999999999999,-5.80000000000001,]

HansenWarren2015_A_strain = [32,32,81,81,118,118,258,258,286,286,337,337,386,386,525,525,]
HansenWarren2015_A_angle  = [48,51,44,35,35,40,1,15,25,25,28,39,1,8,4,11,]

c_red    = sfplt.c_red
c_blue   = sfplt.c_blue
c_green  = '#4daf4a'
c_purple = '#984ea3'
c_orange = '#ff7f00'
c_yellow = '#ffff33'
c_brown  = '#a65628'
c_pink   = '#f781bf'
c_gray   = '#999999'

Skemer_data = {\
    'HansenWarren2015_A':       ('Hansen and Warren (2015)', 's', c_brown, HansenWarren2015_A_strain, HansenWarren2015_A_angle), \
    'Warren2008_A':             ('Warren et al. (2008)', 's', c_green, Warren2008_A_strain, Warren2008_A_angle), \
    'Webber2010_A':             ('Webber et al. (2010)', 's', c_orange, Webber2010_A_strain, Webber2010_A_angle), \
    'Skemer2011_A':             ('Skemer et al. (2011)', 's', c_red, Skemer2011_A_strain, Skemer2011_A_angle), \
    'Hansen2014_A':             ('Hansen et al. (2014)', 's', c_pink, Hansen2014_A_strain, Hansen2014_A_angle), \
    'ZhangKarato1995_1200C_A':  (r'Zhang and Karato (1995), 1200$^\circ$C', 's', c_gray,  ZhangKarato1995_1200C_A_strain, ZhangKarato1995_1200C_A_angle), \
    'ZhangKarato1995_1300C_A':  (r'Zhang and Karato (1995), 1300$^\circ$C', 's', c_purple, ZhangKarato1995_1300C_A_strain, ZhangKarato1995_1300C_A_angle), \
    'Jung2006_A':               ('Jung et al. (2006)', 's', c_blue, Jung2006_A_strain, Jung2006_A_angle), \
}

# experiments with initially isotropic fabric (plot as non-filled markers)
exprdeform = ['ZhangKarato1995_1200C_A', 'ZhangKarato1995_1300C_A', 'Hansen2014_A', 'Skemer2011_A', 'Jung2006_A']

#--------------------
# Plot
#--------------------

### Setup figure

size = 1.6
fig = plt.figure(figsize=(3.1*size,3.5*size))
gs = fig.add_gridspec(ncols=1, nrows=2, height_ratios=[2,1.25])
plt.subplots_adjust(left=0.13, right=1-0.01, top=0.985, bottom=0.370, hspace=0.13)

ax_J = fig.add_subplot(gs[1,0])
ax_phi = fig.add_subplot(gs[0,0], sharex=ax_J)
    
#ax_phi.set_facecolor('#f7f6ee')
#ax_J.set_facecolor('#f7f6ee')
    
sfplt.panellabel(ax_phi, 2, r'\textit{(a)}', frameon=False, fontsize=FS+1.5, bbox=(-0.175,1.09))
sfplt.panellabel(ax_J,   2, r'\textit{(b)}', frameon=False, fontsize=FS+1.5, bbox=(-0.175,1.15))
    
lw = 1.75

color_b = sfplt.c_red
color_n = sfplt.c_blue 

xlims = [0,9]
    
### Observational data

for key in Skemer_data.keys():
    strain = np.divide(Skemer_data[key][-2], 100) # given in percent
    angle  = Skemer_data[key][-1]
    c = Skemer_data[key][2]
    facecolors = 'none' if (key in exprdeform) else c
    ax_phi.scatter(strain, angle, s=6**2, edgecolors=c, facecolors=facecolors, marker=Skemer_data[key][1], label=Skemer_data[key][0], clip_on=1)

### FSE model

dfse = sfdsc.DFSE(dim=3)
dt = sf_time_Lh[1] - sf_time_Lh[0]
ugrad = sf.simpleshear_ugrad(mod['plane'], mod['T'])

FSE_angle = np.zeros(Nt+1)
FSE_angle[0] = 45
for tt in np.arange(1,Nt+1):
    dfse.evolve(ugrad, dt)
#    eigvecs, eigvals, C = dfse.eigenframe()
    C = sfcom.F2C(dfse.F)
    colat = get_m_angles(-C)[0] # colat (index 0) of principal direction (m) (take -C to get antipodal m)
    FSE_angle[tt] = np.rad2deg(sfdsc.colat2lat(colat)) 

ax_phi.plot(f_gammaxz(sf_F_Lh), FSE_angle, '-', color='k', zorder=10, lw=lw, label=r'FSE')
    
#### D-Rex model

I0 = 1 # plot only from this integratation timestep and onwards (shear angle not well defined before this time)
ax_phi.plot(f_gammaxz(drex_F)[I0:], f_shearang(drex_nlm)[I0:], ':', color=color_b, lw=lw, label=r'D-Rex ($b$)')

# fake entry for legend
ax_phi.plot([-1,-1], [0]*2, ':', c=color_n, lw=lw, label=r'D-Rex ($n$)')

### specfab

I0 = 1 # plot only from this integratation timestep and onwards (shear angle not well defined before this time)
ax_phi.plot(f_gammaxz(sf_F_Ll)[I0:],  f_shearang(sf_nlm_Ll)[I0:],  '-',  c=c_orange,      lw=lw, label=r'SDM ($b,n$), $L=%i$'%(Ll))
ax_phi.plot(f_gammaxz(sf_F_Lh)[I0:],  f_shearang(sf_nlm_Lh)[I0:],  '-',  c=color_b,   lw=lw, label=r'SDM ($b,n$), $L=%i$'%(Lh)) 
ax_phi.plot(f_gammaxz(sfd_F_Lh)[I0:], f_shearang(sfd_nlm_Ll)[I0:], '--', c=color_b,   lw=lw, label=r'SDM + $\mathcal{D}$ ($b,n$), $L=%i$'%(Lh)) 

### Labels, legends, ticks, etc.

ax_phi.set_ylabel(r'$\varphi$ ($^\circ$)')

legkwargs = {'handlelength':1.4, 'framealpha':0, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'columnspacing':1.5, 'edgecolor':'k', \
            'borderaxespad':0, 'fontsize':FS-1.7}
handles,labels = ax_phi.get_legend_handles_labels()
hleg = ax_phi.legend(handles,labels, ncol=2, loc='lower right', bbox_to_anchor=(0.92,0.01), bbox_transform=fig.transFigure, **legkwargs)
hleg.get_frame().set_linewidth(0.8)

ylims = [-10,65]
ax_phi.set_yticks(np.arange(0,80,20))
ax_phi.set_yticks(np.arange(-10,80,10), minor=True)
ax_phi.set_xlim(xlims)
ax_phi.set_ylim(ylims)

### ODF insets

if 1:

    geo, prj = sfplt.getprojection(rotation=45, inclination=50)

    def plot_ODF_inset(x, y, nlm, lm, xyarr, title=r'', arr=(0,1)):
        
        W = 0.11 # ax width
        axpos = [x,y, W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        cmap = cmr.get_sub_cmap('Reds', 0, 1) # don't include pure white.
        lvlset = [np.linspace(0.1, 0.6, 6), lambda x,p:'%.1f'%x]
        sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False, nchunk=0)
        sfplt.plotcoordaxes(axin, geo, axislabels='vuxi', color='k', fontsize=FSLEG)
        axin.set_title(title, fontsize=FSLEG)
        
        # Arrow to ODF state
        x_, y_ = xyarr
        ax_phi.annotate("", xy=(x_, y_), xycoords='data', \
                        xytext=(x_+arr[0]/norm, y_+arr[1]/norm), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", facecolor='black'),zorder=20)            

    # ODF inserts

    sf_nlm_trunc = sf_nlm_Lh[:,:sf.L6len] # L=6 truncation for comparrison with D-Rex (only l<=6 coefs derived from D-Rex runs)
    sf_shearang = f_shearang(sf_nlm_trunc)
    drex_shearang = f_shearang(drex_nlm)
    
    I_sf   = lambda Fxz: np.argmin(np.abs(sf_F_Lh[:,0,2]-Fxz))
    I_drex = lambda Fxz: np.argmin(np.abs(drex_F[:,0,2]-Fxz))

    FxzHigh = 5
    FxzLow  = 1

    y0 = 0.43 + 0.26
    
    I = I_sf(FxzLow)
    xyarr = (sf_F_Ll[I,0,2], sf_shearang[I])
    plot_ODF_inset(0.31,0.12+y0, sf_nlm_trunc[I,:], lm_Ll, xyarr, title=r'SDM $(b)$', arr=np.array((0.3,1.5))*0.83)    

    I = I_sf(FxzHigh)
    xyarr = (sf_F_Ll[I,0,2], sf_shearang[I])
    plot_ODF_inset(0.505,0.096+y0, sf_nlm_trunc[I,:], lm_Ll, xyarr, title=r'SDM $(b)$', arr=np.array((-0.125,3.5))*1.05)

    I = I_drex(FxzHigh)
    xyarr = (drex_F[I,0,2], drex_shearang[I])
    plot_ODF_inset(0.67,0.065+y0, drex_nlm[I,:], lm_Ll, xyarr, title=r'D-Rex $(b)$', arr=np.array((+0.3,+4.5))*0.79)

### J-index panel

ax_J.plot(f_gammaxz(sf_F_Lh),  f_J(sf_nlm_Lh),  '-',  c=color_b, lw=lw, label=r'SDM ($b,n$)')
ax_J.plot(f_gammaxz(sfd_F_Lh), f_J(sfd_nlm_Lh), '--', c=color_b, lw=lw, label=r'SDM + $\mathcal{D}$ ($b,n$)')

drex_nlm2, drex_F2 = load_drex_run('simpleShear', 2) # n distr
ax_J.plot(f_gammaxz(drex_F),  f_J(drex_nlm),  ':', c=color_b, lw=lw, label=r'D-Rex ($b$)')
ax_J.plot(f_gammaxz(drex_F2), f_J(drex_nlm2), ':', c=color_n, lw=lw, label=r'D-Rex ($n$)')

# debug
if 0:
    hleg = ax_J.legend(ncol=1, loc=2, **legkwargs)
    hleg.get_frame().set_linewidth(0.8)

ax_J.set_xlabel(r'$\gamma_{xz}$')
plt.setp(ax_phi.get_xticklabels(), visible=False) 

ax_J.set_ylabel(r'$J$')

ylims = [1,9]
ax_J.set_yticks(np.arange(np.round(ylims[0]),ylims[1]+1e-5,2))
ax_J.set_yticks(np.arange(np.round(ylims[0]),ylims[1]+1e-5,1), minor=True)
ax_J.set_ylim(ylims)

ax_J.set_xticks(np.arange(0,12,1))
#ax_J.set_xticks(np.arange(0,12,1), minor=True)
ax_J.set_xlim(xlims)

### Parcel deformation example(s)

ax_J.plot([1]*2, ylims, ':', c='0', lw=1.1, zorder=0)

W = 0.11
ax1 = plt.axes([0.17, 0.460, W,W], projection='3d', zorder=2)
#ax1.patch.set_visible(False)

kwargs = dict(azim=35, axscale=1, axislabels=True, drawinit=True, fonttex=True, fontsize=FSLEG, lw=0.75)
F = sf_F_Lh[np.argmin(np.abs(sf_F_Lh[:,0,2] - 1)),:,:]
sfplt.plotparcel(ax1, F, **kwargs)

### Save fig

fout = 'shear-angle-validation.pdf'
print('Saving %s'%(fout))
plt.savefig(fout, dpi=300)

