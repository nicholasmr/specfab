#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2019-2022

import numpy as np
from numpy import cos, sin, rad2deg, deg2rad
import sys, os, code # code.interact(local=locals())
from scipy.interpolate import interp1d
import scipy.special as sp

sys.path.insert(0, '..')
from specfabpy import specfabpy as sf 

import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.cm as mpl_cm
from matplotlib import colors, ticker, cm 
from matplotlib import rcParams, rc
from matplotlib.offsetbox import AnchoredText
import matplotlib as mpl

lwhatch = 0.9
mpl.rcParams['hatch.linewidth'] = lwhatch
FS = 8.5 + 3.5 + 3.0

rc('font',**{'family':'serif','sans-serif':['Times'],'size':FS})
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage{txfonts} \usepackage{siunitx} \DeclareSIUnit\year{a}'

legkwargs = {'handlelength':1.1, 'framealpha':1.0,'frameon':True, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

def writeSubplotLabel(ax,loc,txt,frameon=True, alpha=1.0, fontsize=FS, pad=0.005, ma='none', bbox=None, zorder=None):
    at = AnchoredText(txt, loc=loc, prop=dict(size=fontsize), frameon=frameon, bbox_to_anchor=bbox, bbox_transform=ax.transAxes)
    at.patch.set_linewidth(0.7)
    if zorder is not None: at.set_zorder(zorder)
    ax.add_artist(at)

#----------------------
# Settings
#----------------------

PRODUCTION = 1

nprime_list = [1,3]

#----------

m, t = np.array([0,0,1]), np.array([1,0,0])
p, q = (m+t)/np.sqrt(2), (m-t)/np.sqrt(2)
mm, mt, pq = np.tensordot(m,m, axes=0), np.tensordot(m,t, axes=0), np.tensordot(p,q, axes=0)

tau0 = 1 # results are independant of the stress magnitude
tau_ps_mm = tau0*(np.identity(3)/3 - mm) 
tau_ss_mt = tau0*(mt + mt.T) 
tau_ss_pq = tau0*(pq + pq.T)

L = 10
lm, nlm_len = sf.init(L) # nlm_len is the number of fabric expansion coefficients (degrees of freedom).
nlm = np.zeros((nlm_len), dtype=np.complex128)
Lrange = np.arange(0,L+1,2) # 0, 2, 4, 6, ...
llrange = [int(ll*(ll+1)/2) for ll in Lrange] # indices 0, 3, 10, ....
for ii, lli in enumerate(llrange):
    nlm[lli] = sp.sph_harm(0, Lrange[ii], 0,0) # expansion coefficients for a delta function

#-----------------------
# MAP
#-----------------------

f = 8 if PRODUCTION else 2
Eca_list = np.logspace(-0,4.1,f*10) # shear along basal plane
Ecc_list = np.logspace(-2,2,f*10) # against basal plane
size = (len(nprime_list), len(Eca_list),len(Ecc_list))
Emm, Emt, Epq = np.zeros(size), np.zeros(size), np.zeros(size)

#-------------

for nn, nprime in enumerate(nprime_list):

    for ii,Eca in enumerate(Eca_list):
        for jj,Ecc in enumerate(Ecc_list):
            Emm[nn,ii,jj] = sf.Evw(nlm, mm,tau_ps_mm, Ecc,Eca,0,nprime) # Evw(a2,a4,a6,a8, vw,tau, Ecc,Eca,alpha,nprime)
            Emt[nn,ii,jj] = sf.Evw(nlm, mt,tau_ss_mt, Ecc,Eca,0,nprime)
            Epq[nn,ii,jj] = sf.Evw(nlm, pq,tau_ss_pq, Ecc,Eca,0,nprime)

        
X = np.array([[ Ecc for Ecc in Ecc_list]   for Eca in Eca_list])
Y = np.array([[ Eca for Ecc in Ecc_list]   for Eca in Eca_list])

#-------------------------    

panelstrs = [r"(a) $n'=1$", r"(b) $n'=3$"]
xysim = [(1,1e4), (1,1e4)]
    
scale = 1.35
plt.figure(figsize=(7.8*scale,2.6*scale))
gs = gridspec.GridSpec(1,2)
gs.update(left=0.07, right=1-0.00, top=0.95, bottom=0.17, wspace=0.2)
ax_list = [plt.subplot(gs[0, 0]),plt.subplot(gs[0, 1])]
cmap = mpl_cm.get_cmap('Blues')
    
for ii,nprime in enumerate(nprime_list):

    print('nprime=%i'%(nprime))

    Emm_map, Emt_map, Epq_map = Emm[ii,:,:], Emt[ii,:,:], Epq[ii,:,:]
    Zpq = np.ma.array(np.divide(Emt_map,Epq_map))
    
    dlvl = 2
    Emt_lvls = np.arange(1,5+1,1)
            
    #--------------------

    ax = ax_list[ii]
    ax.set_xscale('log') 
    ax.set_yscale('log')

    hmap = ax.contourf(X,Y,Emt_map, Emt_lvls, cmap=cmap, vmin=Emt_lvls[0], vmax=Emt_lvls[-1]+0.5, extend='both')
        
    #--------------------
        
    if 1:
        ax1=ax
        
        fmt = ticker.LogFormatterMathtext()
        fmt.create_dummy_axis()
        
        Nneg = 8; N=8
        lvls = [ np.power(10.,n) for n in np.arange(-Nneg,N+1,2) ]
        CS1  = ax1.contour(X,Y,Zpq, lvls, colors='k', linewidths=1.3)
        CS1.collections[0].set_label('$E_{mt}/E_{pq}$')
        
        Nneg = 8; N=8
        lvls = [ np.power(10.,n) for n in np.arange(-Nneg,N+1,2) ]
        CS2  = ax1.contour(X,Y,Emm_map, lvls, linestyles='--',colors='k', linewidths=1.3)
        CS2.collections[0].set_label('$E_{mm}$')

        #----

        def getlblpos(CS,logmid):
        
            label_pos = []
            for line in CS.collections:
                for path in line.get_paths():
                    logvert = np.log10(path.vertices)

                    # find closest point
                    logdist = np.linalg.norm(logvert-logmid, ord=2, axis=1)
                    min_ind = np.argmin(logdist)
                    label_pos.append(10**logvert[min_ind,:])
            return label_pos
            
            
        xmin,xmax,ymin,ymax = Ecc_list[0],Ecc_list[-1], Eca_list[0],Eca_list[-1]

        # work with logarithms for loglog scale
        # middle of the figure:
        logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.5, (np.log10(ymin)+np.log10(ymax))*0.6
        logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.25, (np.log10(ymin)+np.log10(ymax))*0.35
        #
        label_pos_1 = getlblpos(CS1,logmid_1)
        label_pos_2 = getlblpos(CS2,logmid_2)

        # draw labels, hope for the best
        ax1.clabel(CS1, CS1.levels, fmt=fmt, inline_spacing=10, manual=label_pos_1)
        ax1.clabel(CS2, CS2.levels, fmt=fmt, inline_spacing=10, manual=label_pos_2)

    #--------------------

    hsim, = ax.plot(xysim[ii][0],xysim[ii][1], 'X', markeredgecolor='k',markerfacecolor='w', markersize=10)

    legstrs = ['$E_{mt}/E_{pq}$','$E_{mm}$', r'{\fontsize{%i}{%i}\selectfont Simulated}'%(FS-1,FS-1)]
    hlist = [CS1.legend_elements()[0][0], CS2.legend_elements()[0][0], hsim]
    leg=ax.legend(hlist, legstrs, loc=1, bbox_to_anchor=(1,1),fontsize=FS-0.0, **legkwargs); 
    leg.get_frame().set_linewidth(0.7);

    writeSubplotLabel(ax,2, panelstrs[ii],frameon=True, alpha=0.0, fontsize=FS-0.5, pad=-0.05)

    hcb=plt.colorbar(hmap, ax=ax, orientation='vertical', pad=0.06)
    hcb.set_label('$E_{mt}$')
    ax.set_xlabel('$E_{cc}\'$')
    ax.set_ylabel('$E_{ca}\'$')
    ax.set_ylim(Eca_list[[0,-1]])
    
    locmaj = mpl.ticker.LogLocator(base=10,numticks=8) 
    ax.xaxis.set_major_locator(locmaj)
    locmin = mpl.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=8)
    ax.xaxis.set_minor_locator(locmin)
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

plt.savefig('enhancement-factor-nlin-sachs.png', dpi=300)
plt.close()

