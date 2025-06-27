#!/usr/bin/python3
# N. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Fancy version of steadyCPO() class with plotting routines etc. for Rathmann et al. (2025)
"""

import numpy as np
from specfabpy.fenics.steadyCPO import steadyCPO

import matplotlib.pyplot as plt
#import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from matplotlib.offsetbox import AnchoredText
import cmcrameri.cm as cmc

from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=11)

class steadyCPOfancy(steadyCPO):

    kw_vel  = dict(kw_tcf=dict(cmap=cmc.lipari, extend='both', levels=np.logspace(0.5, 3.5, 13), norm='log'), kw_cb=dict(label=r'$u$ (\SI{}{\metre\per\year})'))
    kw_epsE = dict(kw_tcf=dict(cmap=cmc.turku_r, extend='max', levels=np.arange(0, 50+.01, 5)), kw_cb=dict(label=r'$\dot{\epsilon}_{\mathrm{e,iso}}$ (\SI{e-3}{\per\year})'))
    kw_dlam = dict(kw_tcf=dict(cmap=cmc.lapaz_r, extend='max', levels=np.arange(0, 0.8+.01, 0.1)), kw_cb=dict(label=r'$\Delta\lambda$', ticks=np.arange(0, 0.8+.01, 0.2)))
    kw_E    = dict(kw_tcf=dict(cmap='PuOr_r', extend='max', levels=np.arange(0, 4+.01, 0.25), norm='center', vcenter=1), kw_cb=dict(label=r'$E$', ticks=np.arange(0, 4+.01, 1)))

    figsize     = (4,4)
    kw_gs       = dict()
    kw_savefig  = dict(transparent=True, pad_inches=0.05, bbox_inches='tight')
    kw_cax      = dict()
    kw_MODF     = dict(axsize=0.1, xy=[], axloc=[], lbl_bbox=[])
    xticks_major = None
    xticks_minor = None
    yticks_major = None
    yticks_minor = None
    
    kw_leg    = dict(bbox_to_anchor=(-0.07, 1.17), handletextpad=0.5, columnspacing=0.8, handlelength=1.6,)
    color_bcs = ['0.3', 'deeppink'] # isotropic, free
    ls_bcs    = ['-', '--']
    label_bcs = ['Isotropic', 'Free']
    
    color_smax = 'gold'
    label_smax = r'$\vu{z}$ SMAX'
    
    def __init__(self, domain, *args, **kwargs):
        super().__init__(domain)
#        for key, val in kwargs.items(): setattr(self, key, val)

    def newfig(self, **kwargs):
        fig = plt.figure(figsize=self.figsize)
        gs = gridspec.GridSpec(1,2)
        gs.update(**self.kw_gs)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1], sharey=ax1)
        kw_newfig = dict(xticks_major=self.xticks_major, xticks_minor=self.xticks_minor, yticks_major=self.yticks_major, yticks_minor=self.yticks_minor)
        self.setupaxis(ax1, **kw_newfig, **kwargs)
        self.setupaxis(ax2, showyaxis=False, **kw_newfig, **kwargs)
        return (ax1,ax2)

    def savefig(self, num):
        fout = '%s-%i.pdf'%(self.exname, num)
        print('*** Saving %s'%(fout))
        plt.savefig(fout, **self.kw_savefig)

    def plot_inputs(self):
        self.set_inputs()
        ax1, ax2 = self.newfig(mesh=False, boundaries=False, floating=False)
        self.plot_boundaries(ax1)
        self.plot_boundaries(ax2, hidelegend=True)
        self.plot_velocities(ax1, kw_cb=self.kw_vel['kw_cb'], kw_tcf=self.kw_vel['kw_tcf'], kw_cax=self.kw_cax)
        self.plot_strainratemag(ax2, kw_cb=self.kw_epsE['kw_cb'], kw_tcf=self.kw_epsE['kw_tcf'], kw_cax=self.kw_cax)
        self.savefig(0)
      
    def plot_results(self, problem, numerics): 
        self.set_inputs()
        self.set_solution(problem['name'])
        ax1, ax2 = self.newfig(mesh=False, boundaries=False, floating=False)
        if problem['name'] == 'altbc':
            self.color_bcs[0] = self.color_smax
            self.label_bcs[0] = self.label_smax
        self.plot_boundaries(ax1)
        self.plot_boundaries(ax2, hidelegend=True)
        self.plot_E_CAFFE(ax1, kw_cax=self.kw_cax, **self.kw_E)
        self.plot_dlam(   ax2, kw_cax=self.kw_cax, **self.kw_dlam)
        fab = self.femsolution(problem['name'], numerics['L'])
        self.plot_CPOs(ax1, fab)
        self.plot_CPOs(ax2, fab)
        if problem['name'] == 'LROT':       num = 1
        if problem['name'] == 'LROT+DDRX':  num = 2
        if problem['name'] == 'altbc':      num = 3
        self.savefig(num)
        
    def plot_biases(self, problem): 
        self.set_inputs()
        self.set_solution(problem['name'])
        self.c_floating = 'k'
        ax1, ax2 = self.newfig(mesh=False, boundaries=False, floating=False)
        self.plot_boundaries(ax1)
        self.plot_boundaries(ax2, hidelegend=True)
        self.plot_shearfrac(ax1, kw_cax=self.kw_cax)
        self.plot_chi(      ax2, kw_cax=self.kw_cax)
        self.savefig(4 if problem['name'] == 'LROT' else 5)
        
       
    def plot_CPOs(self, ax, fab, onlymarkers=False, lvlmax=0.4, ROTATE_TO_XY=True, dfs=3, cnum='w'):
        geo, prj = sfplt.getprojection(rotation=-90, inclination=(0 if ROTATE_TO_XY else 90))
        for ii, (x,y) in enumerate(self.kw_MODF['xy']):
            if not onlymarkers:
                W = H = self.kw_MODF['axsize'] # ax width
                axin = plt.axes([self.kw_MODF['axloc'][ii][0]-W/2, self.kw_MODF['axloc'][ii][1]-H/2, W, H], projection=prj) 
                axin.set_global()
                nlm = fab.get_nlm(x/self.mapscale, y/self.mapscale)/np.sqrt(4*np.pi) 
                if ROTATE_TO_XY: nlm = fab.sf.rotate_nlm_xz2xy(nlm) # rotate to usual x--y map plane view for SSA models
                sfplt.plotODF(nlm, fab.lm, axin, lvlset=(np.linspace(0.0, lvlmax, 8), lambda x,p:'%.1f'%x), showcb=False, nchunk=None)
                sfplt.plotcoordaxes(axin, geo, axislabels='vuxi', negaxes=False, color='k', fontsize=FS+1.5)
                
            # Markers on map
            num = "%i"%(1+ii)
            peff = [pe.withStroke(linewidth=1.5, foreground='k')]
            if not onlymarkers:
                axin.add_artist(AnchoredText(num, loc=2, prop=dict(color='k', size=FS+1.5), frameon=False, pad=0.3, bbox_to_anchor=self.kw_MODF['lbl_bbox'], bbox_transform=axin.transAxes))
            ax.text(x,y, num, fontsize=FS+dfs, color='w', path_effects=peff, zorder=30, ha='center', va='center') 

    ###########

    def plot_boundaries(self, ax, lw=2, zorder=20, hidelegend=False):
        xb, yb = self.xyboundaries()
        legh, legl = [], []
        for ii in range(2):
            h, = ax.plot(xb[ii]*self.mapscale, yb[ii]*self.mapscale, c=self.color_bcs[ii], lw=lw, ls=self.ls_bcs[ii], zorder=zorder, clip_on=False)
            legh.append(h)
            legl.append(self.label_bcs[ii])
        ax.tricontour(self.triang, self.mask==3, [0.5, 1.5], colors=[self.c_floating,], linewidths=lw, zorder=zorder)
        legh.append(Line2D([0], [0], color=self.c_floating, lw=lw))
        legl.append('Floating')
        if not hidelegend:
            ax.legend(legh, legl, loc='upper left', ncol=3, fancybox=False, frameon=False, **self.kw_leg)

    def xyboundaries(self):
        (coords, bmeshes) = self.bmesh()
        x, y = [], []
        for c in coords:
            xynew = self.reorder_coordinates(c.T)
            x.append(xynew[:,0])
            y.append(xynew[:,1])
        return (x,y)
        
    def reorder_coordinates(self, coords):
        coords = np.array(coords)
        n = len(coords)
        visited = np.zeros(n, dtype=bool)
        path = [0]
        visited[0] = True

        for _ in range(1, n):
            last = path[-1]
            dists = np.linalg.norm(coords - coords[last], axis=1)
            dists[visited] = np.inf
            next_point = np.argmin(dists)
            path.append(next_point)
            visited[next_point] = True

        return self.roll_to_max_distance_endpoints(coords[path])

    def roll_to_max_distance_endpoints(self, coords):
        c = np.vstack((coords, coords[0])).T
        ds = np.sqrt(np.diff(c[0])**2 + np.diff(c[1])**2) # segment lengths
        I = np.argmax(ds)
        return np.roll(coords, -(I+1), axis=0) if I>0 else coords

