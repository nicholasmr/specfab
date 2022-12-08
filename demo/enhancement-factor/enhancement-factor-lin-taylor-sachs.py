#!/usr/bin/python3
# N. M. Rathmann <rathmann@nbi.ku.dk>, 2019-2022

import numpy as np
import sys, os, code # code.interact(local=locals())

sys.path.insert(0, '..')
#from header import *
from specfabpy import specfabpy as sf 

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

# Panels to plot
types = ['Sachs','Taylor','Mixed'] 
#types = ['Taylor']

# High res images?
PRODUCTION = 1

#-------------

nprime = 1 # Linear Sachs
#nprime = 3 # nonlinear Sachs

nprime__Sachs  = nprime
nprime__Taylor = 1

#-------------

m, t = np.array([0,0,1]), np.array([1,0,0])
p, q = (m+t)/np.sqrt(2), (m-t)/np.sqrt(2)
mm, mt, pq = np.tensordot(m,m, axes=0), np.tensordot(m,t, axes=0), np.tensordot(p,q, axes=0)

tau_ps_mm = 1*(np.identity(3)-3*mm) 
tau_ss_mt = 1*(mt + np.transpose(mt)) 
tau_ss_pq = 1*(pq + np.transpose(pq))

L = 4
lm, nlm_len = sf.init(L) # nlm_len is the number of fabric expansion coefficients (degrees of freedom).
a2 = np.tensordot(m,m, axes=0)
a4 = np.tensordot(a2,a2, axes=0)
nlm = sf.a4_to_nlm(a4)

#-----------------------
# Enhancement factor maps
#-----------------------

f = 10 if PRODUCTION else 2

Eca_list = np.logspace(-0,4,f*10) # shear along basal plane
Ecc_list = np.logspace(-2,2,f*10) # against basal plane
size = (len(Eca_list),len(Ecc_list))

if nprime__Sachs == 3:  alpha_list = np.logspace(-3.4,0,f*10) # Sachs--Taylor weight 
else:                   alpha_list = np.logspace(-3,0,f*10)   # Sachs--Taylor weight 

Emm_Sachs,  Emt_Sachs,  Epq_Sachs  = np.zeros(size), np.zeros(size), np.zeros(size)
Emm_Taylor, Emt_Taylor, Epq_Taylor = np.zeros(size), np.zeros(size), np.zeros(size)
Emm_alpha,  Emt_alpha,  Epq_alpha  = np.zeros(size), np.zeros(size), np.zeros(size)

#-------------

#print(psi.Evw(mt,tau_ss_mt, 3,[1],[1], moments))
#print()
#print('max E_mt (n=1): %f'%(psi.Evw(mt,tau_ss_mt, 1,[0.1],[1e8], moments)[0][0]))
#print('max E_mt (n=3): %f'%(psi.Evw(mt,tau_ss_mt, 3,[0.1],[1e8], moments)[0][0]))

#-------------

for ii,Eca in enumerate(Eca_list):

    for jj,Ecc in enumerate(Ecc_list):
        Emm_Sachs[ii,jj] = sf.Evw(nlm, mm,tau_ps_mm, Ecc,Eca,0,nprime)
        Emt_Sachs[ii,jj] = sf.Evw(nlm, mt,tau_ss_mt, Ecc,Eca,0,nprime)
        Epq_Sachs[ii,jj] = sf.Evw(nlm, pq,tau_ss_pq, Ecc,Eca,0,nprime)
        #       
        Emm_Taylor[ii,jj] = sf.Evw(nlm, mm,tau_ps_mm, Ecc,Eca,1,nprime)
        Emt_Taylor[ii,jj] = sf.Evw(nlm, mt,tau_ss_mt, Ecc,Eca,1,nprime)
        Epq_Taylor[ii,jj] = sf.Evw(nlm, pq,tau_ss_pq, Ecc,Eca,1,nprime)
        
    for kk,alpha in enumerate(alpha_list):
        Emm_alpha[ii,kk] = sf.Evw(nlm, mm,tau_ps_mm, 1,Eca,alpha,nprime)
        Emt_alpha[ii,kk] = sf.Evw(nlm, mt,tau_ss_mt, 1,Eca,alpha,nprime)
        Epq_alpha[ii,kk] = sf.Evw(nlm, pq,tau_ss_pq, 1,Eca,alpha,nprime)
        
Xe = np.array([[ Ecc for Ecc in Ecc_list]   for Eca in Eca_list])
Xa = np.array([[ alp for alp in alpha_list] for Eca in Eca_list])
Y  = np.array([[ Eca for Ecc in Ecc_list]   for Eca in Eca_list])

#-----------------------
# Plot maps
#-----------------------

panelstrs = [r'{\LARGE \bf a}\,\, $\alpha=0$', r'{\LARGE \bf b}\,\, $\alpha=1$', r'{\LARGE \bf c}\,\, $E_{cc} = 1$']

for ii,TYPE in enumerate(types):

    print(TYPE)

    if TYPE is 'Sachs':
        nprime = nprime__Sachs;
        Emm_map, Emt_map, Epq_map = Emm_Sachs, Emt_Sachs, Epq_Sachs
        X = Xe; xlbl='$E_{cc}\'$';
        cmap = mpl_cm.get_cmap('Blues')
        if nprime__Sachs == 3:  
            dlvl = 2
            Emt_lvls = np.arange(1,5+1,1)
        else:                   
            dlvl = 1
            Emt_lvls = np.arange(1,2.5+0.5,0.5)
            
    if TYPE is 'Taylor':
        nprime = nprime__Taylor;
        Emm_map, Emt_map, Epq_map = Emm_Taylor, Emt_Taylor, Epq_Taylor
        X = Xe; xlbl='$E_{cc}\'$';
        cmap = mpl_cm.get_cmap('Greens')
        dlvl = 1
        Emt_lvls = np.arange(-0,4+1,1)
        
    if TYPE is 'Mixed':
        n = None;
        Emm_map, Emt_map, Epq_map = Emm_alpha, Emt_alpha, Epq_alpha
        X = Xa; xlbl=r'$\alpha$';
        cmap = mpl_cm.get_cmap('Oranges')
        dlvl = 1
        Emt_lvls = np.arange(-0,3+1,1)

    Zpq = np.ma.array(np.divide(Emt_map,Epq_map))
    
    #-------------

    scale=1.35
    plt.figure(figsize=(2.95*scale,3.6*scale))
    gs = gridspec.GridSpec(1,1)
    a = 0.18
    gs.update(left=a, right=1-a*0.3, top=0.95, bottom=0.16, wspace=0.13)
    ax1 = plt.subplot(gs[0, 0])
    ax1.set_xscale('log') 
    ax1.set_yscale('log')

    if TYPE is 'Sachs':
        hmap = ax1.contourf(X,Y,Emt_map, Emt_lvls, cmap=cmap, vmin=Emt_lvls[0], vmax=Emt_lvls[-1]+0.5, extend='both')
    else:
        hmap = ax1.contourf(X,Y,Emt_map, np.power(10., Emt_lvls), cmap=cmap, norm=colors.LogNorm(), extend='both') # https://matplotlib.org/3.3.2/gallery/images_contours_and_fields/contourf_log.html
        
    if 1:
    
        fmt = ticker.LogFormatterMathtext()
        fmt.create_dummy_axis()
        
        Nneg = 8; N=8
        if TYPE is "Mixed": Nneg=-1
    #    lvls = np.logspace(-Nneg,N,N+Nneg+1)
        lvls = [ np.power(10.,n) for n in np.arange(-Nneg,N+1,dlvl) ]
        CS1  = ax1.contour(X,Y,Zpq, lvls, colors='k', linewidths=1.3)
        CS1.collections[0].set_label('$E_{mt}/E_{pq}$')
        
        Nneg = 8; N=8
    #    lvls = np.logspace(-Nneg,N,N+Nneg+1)
        lvls = [ np.power(10.,n) for n in np.arange(-Nneg,N+1,dlvl) ]
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

        if TYPE is "Sachs":
            logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.5, (np.log10(ymin)+np.log10(ymax))*0.6
            logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.25, (np.log10(ymin)+np.log10(ymax))*0.35

        if TYPE is "Taylor":
            logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.3, (np.log10(ymin)+np.log10(ymax))*0.3
            logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.25, (np.log10(ymin)+np.log10(ymax))*0.40

        if TYPE is "Mixed":
            if nprime__Sachs == 3:  
                logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.73, (np.log10(ymin)+np.log10(ymax))*0.55
                logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.26, (np.log10(ymin)+np.log10(ymax))*0.35
            else:                   
                logmid_1 = (np.log10(xmin)+np.log10(xmax))*0.85, (np.log10(ymin)+np.log10(ymax))*0.55
                logmid_2 = (np.log10(xmin)+np.log10(xmax))*0.625, (np.log10(ymin)+np.log10(ymax))*0.5

        #
        label_pos_1 = getlblpos(CS1,logmid_1)
        label_pos_2 = getlblpos(CS2,logmid_2)

        # draw labels
        inline_spacing=12
#        inline_spacing=20
        ax1.clabel(CS1, CS1.levels, fmt=fmt, inline_spacing=inline_spacing, manual=label_pos_1)
        ax1.clabel(CS2, CS2.levels, fmt=fmt, inline_spacing=inline_spacing, manual=label_pos_2)

    #-------------

    if TYPE is "Mixed":
        xysim = (1.25e-2,1e3)
        hsim, = ax1.plot(xysim[0],xysim[1], 'X', markeredgecolor='k',markerfacecolor='w', markersize=10)

    hlist = [CS1.legend_elements()[0][0], CS2.legend_elements()[0][0]]
    hlbls = ['$E_{mt}/E_{pq}$','$E_{mm}$']
    if TYPE is "Mixed":
        hlist += [hsim]
        hlbls += [r'{\fontsize{%i}{%i}\selectfont Simulated}'%(FS-1,FS-1)]
    leg=ax1.legend(hlist, hlbls, loc=1, bbox_to_anchor=(1,1),fontsize=FS-0.0, **legkwargs); 
    leg.get_frame().set_linewidth(0.7);

    writeSubplotLabel(ax1,2, panelstrs[ii],frameon=True, alpha=0.0, fontsize=FS-0.5, pad=-0.05)

    hcb=plt.colorbar(hmap, orientation='horizontal', pad=0.18)
    hcb.set_label('$E_{mt}$')
    ax1.set_xlabel(xlbl)
    ax1.set_ylabel('$E_{ca}\'$')
    ax1.set_ylim(Eca_list[[0,-1]])
    
    locmaj = mpl.ticker.LogLocator(base=10,numticks=8) 
    ax1.xaxis.set_major_locator(locmaj)
    locmin = mpl.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=8)
    ax1.xaxis.set_minor_locator(locmin)
    ax1.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    plt.savefig('Evw_%s.pdf'%(TYPE), dpi=200)
    plt.close()
    
if len(types) == 3:
    flist = "Evw_Sachs.pdf Evw_Taylor.pdf Evw_Mixed.pdf"
    fout = "enhancement-factor-lin-taylor-sachs"
    os.system('pdfjam --nup 3x1 --trim "0 0 0em 0" %s --outfile %s.pdf;'%(flist,fout))    
    os.system('pdfcrop --margins 2 %s.pdf %s.pdf;'%(fout,fout))
    os.system('convert -density 200 -quality 100 +profile "icc" -flatten %s.pdf %s.png'%(fout,fout))
    os.system('rm %s'%(flist))
    os.system('rm %s.pdf'%(fout))


