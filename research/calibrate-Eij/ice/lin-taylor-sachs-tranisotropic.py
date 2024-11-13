#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2019-2023

import sys, os, code # code.interact(local=locals())

import numpy as np
from specfabpy import specfab as sf 
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt

sys.path.append('..')
from localheader import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as mpl_cm
from matplotlib import colors, ticker 
import matplotlib as mpl

import warnings
warnings.filterwarnings("ignore")

lwhatch = 0.9
mpl.rcParams['hatch.linewidth'] = lwhatch

FS = 8.5 + 3.5 + 3.0
sfplt.setfont_tex(fontsize=FS)

#----------------------
# Settings
#----------------------

types = ['Sachs','Taylor','Mixed'] 
#types = ['Taylor']

PRODUCTION = 1

#-----------------------
# Generate maps
#-----------------------

f = 10 if PRODUCTION else 2

Eca_list   = np.logspace(-0,4,f*10) # shear along basal plane
Ecc_list   = np.logspace(-2,2,f*10) # against basal plane
alpha_list = np.logspace(-3,0,f*10) # Sachs--Taylor weight 

n_grain = 1

Emm_map_Taylor, Emt_map_Taylor, Epq_map_Taylor, Xe, Y = Eij_maps_tranisotropic(Ecc_list,   Eca_list, 'Taylor', n_grain=n_grain)
Emm_map_Sachs,  Emt_map_Sachs,  Epq_map_Sachs,  Xe, Y = Eij_maps_tranisotropic(Ecc_list,   Eca_list, 'Sachs',  n_grain=n_grain)
Emm_map_Mixed,  Emt_map_Mixed,  Epq_map_Mixed,  Xa, Y = Eij_maps_tranisotropic(alpha_list, Eca_list, 'Mixed',  n_grain=n_grain, Ecc=1)

#-----------------------
# Plot
#-----------------------

panelstrs = [r'\textit{(a)}\, $\alpha=0$', r'\textit{(b)}\, $\alpha=1$', r'\textit{(c)}\, $E_{cc}^\prime = 1$']

for ii,TYPE in enumerate(types):

    print(TYPE)

    ### Determine plot bounds

    if TYPE == 'Sachs':
        Emm_map, Emt_map, Epq_map = Emm_map_Sachs, Emt_map_Sachs, Epq_map_Sachs
        X = Xe
        xlbl = r'$E_{cc}^\prime$'
        cmap = mpl_cm.get_cmap('Blues')
        dlvl = 1
        Emt_lvls = np.arange(1,2.5+0.5,0.5)
            
    if TYPE == 'Taylor':
        Emm_map, Emt_map, Epq_map = Emm_map_Taylor, Emt_map_Taylor, Epq_map_Taylor
        X = Xe
        xlbl = r'$E_{cc}^\prime$'
        cmap = mpl_cm.get_cmap('Greens')
        dlvl = 1
        Emt_lvls = np.arange(-0,4+1,1)
        
    if TYPE == 'Mixed':
        Emm_map, Emt_map, Epq_map = Emm_map_Mixed, Emt_map_Mixed, Epq_map_Mixed
        X = Xa
        xlbl = r'$\alpha$'
        cmap = mpl_cm.get_cmap('Oranges')
        dlvl = 1
        Emt_lvls = np.arange(-0,3+1,1)

    Zpq = np.ma.array(np.divide(Emt_map,Epq_map))
    
    ### Setup figure

    scale=1.35
    plt.figure(figsize=(2.95*scale,3.6*scale))
    gs = gridspec.GridSpec(1,1)
    a = 0.18
    gs.update(left=a, right=1-a*0.3, top=0.95, bottom=0.16, wspace=0.13)
    ax1 = plt.subplot(gs[0, 0])
    ax1.set_xscale('log') 
    ax1.set_yscale('log')

    legkwargs = {'handlelength':1.1, 'framealpha':1.0,'frameon':True, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

    ### Plot

    if TYPE == 'Sachs':
        hmap = ax1.contourf(X,Y,Emt_map, Emt_lvls, cmap=cmap, vmin=Emt_lvls[0], vmax=Emt_lvls[-1]+0.5, extend='both')
    else:
        hmap = ax1.contourf(X,Y,Emt_map, np.power(10., Emt_lvls), cmap=cmap, norm=colors.LogNorm(), extend='both') # https://matplotlib.org/3.3.2/gallery/images_contours_and_fields/contourf_log.html
        
    if 1:
    
        fmt = ticker.LogFormatterMathtext()
        fmt.create_dummy_axis()
        
        Nneg = 8; N=8
        if TYPE == "Mixed": Nneg=-1
        lvls = [ np.power(10.,n) for n in np.arange(-Nneg,N+1,dlvl) ]
        CS1  = ax1.contour(X,Y,Zpq, lvls, colors='k', linewidths=1.3)
        CS1.collections[0].set_label('$E_{mt}/E_{pq}$')
        
        Nneg = 8; N=8
        lvls = [ np.power(10.,n) for n in np.arange(-Nneg,N+1,dlvl) ]
        if 0: # DEBUG: for comparing with orthotropic case 
            lvls = [5e-2,1e-1,2.5e-1,5e-1,7.5e-1] 
            fmt = '%.2f'
        CS2  = ax1.contour(X,Y,Emm_map, lvls, linestyles='--',colors='k', linewidths=1.3)
        CS2.collections[0].set_label('$E_{mm}$')

        #-------------

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
            
        xmin,xmax,ymin,ymax = plt.axis()

        if TYPE == "Sachs":
            logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.5, (np.log10(ymin)+np.log10(ymax))*0.6
            logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.25, (np.log10(ymin)+np.log10(ymax))*0.35

        if TYPE == "Taylor":
            logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.3, (np.log10(ymin)+np.log10(ymax))*0.3
            logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.25, (np.log10(ymin)+np.log10(ymax))*0.40

        if TYPE == "Mixed":
            logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.85, (np.log10(ymin)+np.log10(ymax))*0.55
            logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.625, (np.log10(ymin)+np.log10(ymax))*0.5

        label_pos_1 = getlblpos(CS1,logmid_1)
        label_pos_2 = getlblpos(CS2,logmid_2)

        inline_spacing=12
#        inline_spacing=20
        ax1.clabel(CS1, CS1.levels, fmt=fmt, inline_spacing=inline_spacing, manual=label_pos_1)
        ax1.clabel(CS2, CS2.levels, fmt=fmt, inline_spacing=inline_spacing, manual=label_pos_2)

    #-------------

    if TYPE == "Mixed":
        xysim = (1.25e-2, 1e3)
        hsim, = ax1.plot(xysim[0],xysim[1], 'X', markeredgecolor='k',markerfacecolor='w', markersize=10)

    hlist = [CS1.legend_elements()[0][0], CS2.legend_elements()[0][0]]
    hlbls = ['$E_{mt}/E_{pq}$','$E_{mm}$']
    if TYPE == "Mixed":
        hlist += [hsim]
        hlbls += [r'{\fontsize{%i}{%i}\selectfont Simulated}'%(FS-1,FS-1)]
    leg=ax1.legend(hlist, hlbls, loc=1, bbox_to_anchor=(1,1),fontsize=FS-0.0, **legkwargs); 
    leg.get_frame().set_linewidth(0.7);

    sfplt.panellabel(ax1,2, panelstrs[ii],frameon=True, alpha=0.0, fontsize=FS-0.5, pad=0.4)

    hcb=plt.colorbar(hmap, orientation='horizontal', pad=0.18)
    hcb.set_label('$E_{mt}$')
    ax1.set_xlabel(xlbl)
    ax1.set_ylabel('$E_{ca}\'$')
    ax1.set_ylim(Eca_list[[0,-1]])
    ax1.set_xlim(alpha_list[[0,-1]] if TYPE == 'Mixed' else Ecc_list[[0,-1]])
    
    locmaj = mpl.ticker.LogLocator(base=10,numticks=8) 
    ax1.xaxis.set_major_locator(locmaj)
    locmin = mpl.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=8)
    ax1.xaxis.set_minor_locator(locmin)
    ax1.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    plt.savefig('Evw_%s.pdf'%(TYPE), dpi=200)

if len(types) == 3:
    flist = "Evw_Sachs.pdf Evw_Taylor.pdf Evw_Mixed.pdf"
    fout = "lin-taylor-sachs-tranisotropic"
    os.system('pdfjam --nup 3x1 --trim "0 0 0em 0" %s --outfile %s.pdf;'%(flist,fout))    
    os.system('pdfcrop --margins 2 %s.pdf %s.pdf;'%(fout,fout))
    os.system('convert -density 200 -quality 100 +profile "icc" -flatten %s.pdf %s.png'%(fout,fout))
    os.system('rm %s'%(flist))
#    os.system('rm %s.pdf'%(fout))

