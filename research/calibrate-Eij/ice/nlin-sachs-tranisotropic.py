#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2019-2023

import sys, os, code # code.interact(local=locals())

import numpy as np

from specfabpy import specfab as sf 
from specfabpy import plotting as sfplt

sys.path.append('..')
from localheader import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as mpl_cm
from matplotlib import ticker 
import matplotlib as mpl

import warnings
warnings.filterwarnings("ignore")

lwhatch = 0.9
mpl.rcParams['hatch.linewidth'] = lwhatch

FS = 8.5 + 3.5 + 3.0
sfplt.setfont_tex(fontsize=FS)
legkwargs = {'handlelength':1.1, 'framealpha':1.0,'frameon':True, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

#----------------------
# Settings
#----------------------

PRODUCTION = 1

#-----------------------
# Generate maps
#-----------------------

f = 8 if PRODUCTION else 2
Eca_list = np.logspace(-0,4.1,f*10) # shear along basal plane
Ecc_list = np.logspace(-2,2,f*10) # against basal plane
Emm_map_n1, Emt_map_n1, Epq_map_n1, X, Y = Eij_maps_tranisotropic(Ecc_list, Eca_list, 'Sachs', n_grain=1)
Emm_map_n3, Emt_map_n3, Epq_map_n3, X, Y = Eij_maps_tranisotropic(Ecc_list, Eca_list, 'Sachs', n_grain=3)

#-----------------------
# Plot
#-----------------------    

panelstrs = [r"\textit{(a)}\, $n'=1$", r"\textit{(b)}\, $n'=3$"]
xysim = [(1,1e4), (1,1e4)]
cmap = mpl_cm.get_cmap('Blues')
    
scale = 1.35
plt.figure(figsize=(7.8*scale,2.6*scale))
gs = gridspec.GridSpec(1,2)
gs.update(left=0.07, right=1-0.00, top=0.95, bottom=0.17, wspace=0.2)
ax_list = [plt.subplot(gs[0, 0]),plt.subplot(gs[0, 1])]
    
for ii,n_grain in enumerate([1,3]):

    print('n_grain=%i'%(n_grain))

    if n_grain==1: Emm_map, Emt_map, Epq_map = Emm_map_n1, Emt_map_n1, Epq_map_n1
    if n_grain==3: Emm_map, Emt_map, Epq_map = Emm_map_n3, Emt_map_n3, Epq_map_n3
    Zpq = np.ma.array(np.divide(Emt_map,Epq_map))
    
    dlvl = 2
    Emt_lvls = np.arange(1,5+1,1)
            
    ax = ax_list[ii]
    ax.set_xscale('log') 
    ax.set_yscale('log')

    hmap = ax.contourf(X,Y,Emt_map, Emt_lvls, cmap=cmap, vmin=Emt_lvls[0], vmax=Emt_lvls[-1]+0.5, extend='both')
        
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

    sfplt.panellabel(ax,2, panelstrs[ii], frameon=True, alpha=0.0, fontsize=FS-0.5, pad=0.4)

    hcb=plt.colorbar(hmap, ax=ax, orientation='vertical', pad=0.06)
    hcb.set_label('$E_{mt}$')
    ax.set_xlabel(r'$E_{cc}^\prime$')
    ax.set_ylabel(r'$E_{ca}^\prime$')
    ax.set_ylim(Eca_list[[0,-1]])
    ax.set_xlim(Ecc_list[[0,-1]])
    
    locmaj = mpl.ticker.LogLocator(base=10,numticks=8) 
    ax.xaxis.set_major_locator(locmaj)
    locmin = mpl.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=8)
    ax.xaxis.set_minor_locator(locmin)
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

plt.savefig('nlin-sachs-tranisotropic.png', dpi=250)

