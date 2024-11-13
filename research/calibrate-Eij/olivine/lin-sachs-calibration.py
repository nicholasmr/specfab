#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

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

PRODUCTION = 1

#-----------------------
# Generate maps
#-----------------------

f = 10 if PRODUCTION else 2

Enb_list   = np.logspace(-0,4,f*10) # shear along basal plane
Eii_list   = np.logspace(-2,2,f*10) # against basal plane
alpha_list = np.logspace(-3,0,f*10) # Sachs--Taylor weight 

n_grain = 1

Emm_map_Sachs,  Emt_map_Sachs,  Epq_map_Sachs,  Xe, Y = Eij_maps_orthotropic(Eii_list,   Enb_list, 'Sachs',  n_grain=n_grain)

#-----------------------
# Plot maps
#-----------------------

xysim = (0.3,1e1) # values used for simulation in mixed (alpha) map

legkwargs = {'handlelength':1.1, 'framealpha':1.0,'frameon':True, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

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

scale=1.2
plt.figure(figsize=(7.2*scale,3.0*scale))
gs = gridspec.GridSpec(1,2, width_ratios=[1.2,1])
a = 0.09
gs.update(left=a, right=1-a/4, top=0.92, bottom=0.17, wspace=0.2)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])

### Sachs map 

ax1.set_xscale('log') 
ax1.set_yscale('log')

Emm_map, Emt_map = Emm_map_Sachs, Emt_map_Sachs
X = Xe
xlbl = '$E_{ii}^\prime$'
cmap = mpl_cm.get_cmap('Blues')
dlvl = 1
Emt_lvls = np.arange(1,5+0.5,1)

hmap = ax1.contourf(X,Y,Emt_map, Emt_lvls, cmap=cmap, vmin=Emt_lvls[0], vmax=Emt_lvls[-1]+0.5, extend='both')

fmt = ticker.LogFormatterMathtext(labelOnlyBase=False)
fmt.create_dummy_axis()

Nneg = 8; N=8
lvls = [ np.power(10.,n) for n in np.arange(-Nneg,N+1,dlvl) ]
CS2  = ax1.contour(X,Y,Emm_map, lvls, linestyles='--', colors='k', linewidths=1.3)

xmin,xmax,ymin,ymax = ax1.axis()
logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.25, (np.log10(ymin)+np.log10(ymax))*0.35
label_pos_2 = getlblpos(CS2,logmid_2)

inline_spacing=12
ax1.clabel(CS2, CS2.levels, fmt=fmt, inline_spacing=inline_spacing, manual=label_pos_2)

hlist = [CS2.legend_elements()[0][0],]
hlbls = ['$E_{11}$',]

leg=ax1.legend(hlist, hlbls, loc=1, bbox_to_anchor=(1,1),fontsize=FS-0.0, **legkwargs); 
leg.get_frame().set_linewidth(0.8);

hcb=plt.colorbar(hmap, orientation='vertical', pad=0.05)
hcb.set_label('$E_{12}$')
ax1.set_xlabel(xlbl)
ax1.set_ylabel('$E_{nb}\'$')
ax1.set_ylim(Enb_list[[0,-1]])

locmaj = mpl.ticker.LogLocator(base=10, numticks=8) 
ax1.xaxis.set_major_locator(locmaj)
locmin = mpl.ticker.LogLocator(base=10.0, subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), numticks=8)
ax1.xaxis.set_minor_locator(locmin)
ax1.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

### Eii = 1 line 

xlims = [Y.min(), Y.max()]
ylims = [0, 5+1e-3]

#    ax2.semilogx(xlims, [5]*2, ':', c=sfplt.c_dblue)
#    ax2.semilogx(xlims, [0.5]*2, ':', c='k')

I = np.argmin(np.abs(Eii_list-1))
ax2.semilogx(Y[:,I], Emt_map[:,I], '-',  c='k', label='$E_{12}$')
ax2.semilogx(Y[:,I], Emm_map[:,I], '--', c='k', label='$E_{11}$')

def plot_Enb_ref(Enb_ref, c='0.5', lw=1.75):
    Inb = np.argmin(np.abs(Enb_list-Enb_ref))
    ax2.semilogx([Enb_ref]*2, ylims, ':', c=c, lw=lw)
    Emt_ref, Emm_ref = Emt_map[Inb,I], Emm_map[Inb,I]
    ax2.plot(Enb_ref, Emt_ref, 'o', ms=9, markeredgecolor='k', markerfacecolor="None", markeredgewidth=1.5)
    ax2.plot(Enb_ref, Emm_ref, 'o', ms=9, markeredgecolor='k', markerfacecolor="None", markeredgewidth=1.5)
    ax2.semilogx(xlims, [Emt_ref]*2, ':', c=c, lw=lw)
    ax2.semilogx(xlims, [Emm_ref]*2, ':', c=c, lw=lw)
    print('Enb=%i: Emt=%.4f'%(Enb_ref, Emt_ref))
    print('Enb=%i: Emm=%.4f'%(Enb_ref, Emm_ref))
        
plot_Enb_ref(100, c='#b35806')
plot_Enb_ref(20, c='#542788')

ax2.set_yticks(np.arange(*ylims,1))
ax2.set_yticks(np.arange(*ylims,0.1), minor=True)

ax2.set_ylabel('$E_{ij}$')
ax2.set_xlabel('$E_{nb}\'$')
ax2.set_xlim(xlims)
ax2.set_ylim(ylims)

locmaj = mpl.ticker.LogLocator(base=10, numticks=8) 
ax2.xaxis.set_major_locator(locmaj)
locmin = mpl.ticker.LogLocator(base=10.0, subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), numticks=8)
ax2.xaxis.set_minor_locator(locmin)
ax2.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

leg = ax2.legend(loc=1, **legkwargs)
leg.get_frame().set_linewidth(0.8);
sfplt.panellabel(ax2, 4, r'$E_{ij}^\prime = 1 \forall i,j\neq n,b$', frameon=False, alpha=0.0, fontsize=FS-1, pad=0.0, bbox=(1,0.05))

sfplt.panellabel(ax1, 2, r'\textit{(a)}', frameon=False, fontsize=FS+2, bbox=(-0.30,1.14))
sfplt.panellabel(ax2, 2, r'\textit{(b)}', frameon=False, fontsize=FS+2, bbox=(-0.24,1.14))
fout = "lin-sachs-orthotropic"
plt.savefig('%s.pdf'%(fout), dpi=200)
plt.close()
#    os.system('pdfcrop %s.pdf %s.pdf'%(fout,fout))
os.system('pdftocairo -singlefile -r 250 -png %s.pdf %s'%(fout,fout))

