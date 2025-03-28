#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Common routines for Ross and PIG plots in Rathmann et al. (2025)
"""

from experiment import *

import numpy as np
from dolfin import *

from specfabpy import specfab as sf
from specfabpy.fenics.ice import IceFabric
from specfabpy import plotting as sfplt

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import matplotlib.patheffects as pe
import cmcrameri.cm as cmc

FS = sfplt.setfont_tex(fontsize=11)

class R25(Experiment):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.veltransform=lambda x:x

    def makefig(self, lbla, lblb, boundaries=True):
        fig = plt.figure(figsize=self.figsize)
        gs = gridspec.GridSpec(1,2)
        gs.update(**self.kwargs_gs)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1], sharey=ax1)

        for ii, ax in enumerate((ax1,ax2)):
            ax.add_patch(plt.Rectangle((self.xlims[0],self.ylims[0]), np.diff(self.xlims)[0], np.diff(self.ylims)[0], color='0.85'))
            ax.axis('square')
            ax.set_xlabel(r'$x$ ($\SI{}{\kilo\metre}$)')
            ax.set_xticks(self.xticks_major)
            if self.xticks_minor is not None: ax.set_xticks(self.xticks_minor, minor=True)
            ax.set_xlim(self.xlims)

            if ii == 0:
                ax.set_ylabel(r'$y$ ($\SI{}{\kilo\metre}$)')
                ax.set_yticks(self.yticks_major)
                if self.yticks_minor is not None: ax.set_yticks(self.yticks_minor, minor=True)
                ax.set_ylim(self.ylims)
            else:
                ax.tick_params('y', labelleft=False)
                ax.set_ylim(self.ylims)
     
        sfplt.panellabel(ax1, 2, r'\fontsize{13}{13}\selectfont\textbf{%s}'%(lbla), **self.kwargs_lbl)
        sfplt.panellabel(ax2, 2, r'\fontsize{13}{13}\selectfont\textbf{%s}'%(lblb), **self.kwargs_lbl)

        return ax1,ax2

    def savefig(self, num):
        fout = '%s/%s-%i.pdf'%(self.exname, self.exname, num)
        print('*** Saving %s'%(fout))
        plt.savefig(fout, transparent=True, dpi=200)
        #os.system('pdfcrop %s %s'%(fout, fout))

    #### Load data

    def plot_velocities(self, lbl1='a', lbl2='b'):

        # FEM gridded data
        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        coords,cells, ux,uy,umag,epsE, S,B,H,mask = self.get_np_inputs()
        triang = self.get_triang(coords, cells, scale=self.lenscale)

        ### Plot observed velocities *as stored on FE mesh*

        ax1,ax2 = self.makefig(lbl1, lbl2)

        ax = ax1
        cs = ax.tricontourf(triang, self.veltransform(self.velscale*umag), levels=self.lvls_vel, extend='both', locator=None, cmap=cmc.lipari)
        hcb = plt.colorbar(cs, **self.kwargs_cb)
        hcb.set_label(r'%s (\SI{}{\metre\per\year})'%(self.clabel_vel))
        self.plot_domain(ax, triang, mask, scale=self.lenscale, legend=True)
#        ax1.triplot(triang, **{'lw':0.075, 'color':'0.5', 'alpha':0.8})
#        self.plot_stream(ax, sc*X2, sc*Y2, UX2, UY2)

        ax = ax2
        F = 1e3*self.velscale * epsE
        cs = ax.tricontourf(triang, F, levels=self.lvls_strainrate, extend='max', locator=None, cmap=cmc.turku_r)
        hcb = plt.colorbar(cs, **self.kwargs_cb)
        hcb.set_label(r'$\dot{\epsilon}_{\mathrm{iso}}$ (\SI{e-3}{\per\year})'%())
        self.plot_domain(ax, triang, mask, scale=self.lenscale)

        self.savefig(0)

    def plot_results(self, FP, lbl1='a', lbl2='b', lvlsE_max=3, kwargs_MODF=dict()):

        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        coords,cells, ux,uy,umag,epsE, S,B,H,mask = self.get_np_inputs()
        coords,cells, mi,lami,E_CAFFE = self.get_np_solution(FP)
        triang = self.get_triang(coords, cells, scale=self.lenscale)
        
        ax1, ax2 = self.makefig(lbl1, lbl2)
        self._plot_E(ax1, triang, E_CAFFE, lvlsE_max=lvlsE_max)
        self.plot_domain(ax1, triang, mask, scale=self.lenscale, legend=True)
        
        self._plot_dlam(ax2, triang, lami)
        self.plot_domain(ax2, triang, mask, scale=self.lenscale)
        
        fab = self.get_FE_solution(FP)
        self.plot_CPOs(ax1, ax2, fab, **kwargs_MODF)

        self.savefig(1 if FP=='LROT' else 2)
        
    def plot_altbc_compare(self, FP, lvlsE_max=3, kwargs_MODF=dict()):
    
        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        coords,cells, ux,uy,umag,epsE, S,B,H,mask = self.get_np_inputs()
        triang = self.get_triang(coords, cells, scale=self.lenscale)
        
        ### Iso BC
        
        coords,cells, mi,lami,E_CAFFE = self.get_np_solution(FP)        
        ax1, ax2 = self.makefig('a', 'b')
        self._plot_E(ax1, triang, E_CAFFE, lvlsE_max=lvlsE_max)
        self.plot_domain(ax1, triang, mask, scale=self.lenscale, legend=True)
        self._plot_dlam(ax2, triang, lami)
        self.plot_domain(ax2, triang, mask, scale=self.lenscale)
        self.plot_CPOs(ax1, ax2, self.get_FE_solution(FP), **kwargs_MODF)
        self.savefig(10)
        
        ### Alt BC
    
        coords,cells, mi,lami,E_CAFFE = self.get_np_solution(FP, suffix='-altbc')
        ax1, ax2 = self.makefig('c', 'd')
        self._plot_E(ax1, triang, E_CAFFE, lvlsE_max=lvlsE_max)
        self.plot_domain(ax1, triang, mask, scale=self.lenscale, legend=True, lbl_inflow=r'$\vu{z}$ SMAX', c_inflow='gold')
        self._plot_dlam(ax2, triang, lami)
        self.plot_domain(ax2, triang, mask, scale=self.lenscale)
        self.plot_CPOs(ax1, ax2, self.get_FE_solution(FP, suffix='-altbc'), **kwargs_MODF)
        self.savefig(11)


    def _plot_E(self, ax, triang, E, lvlsE_max=3):
#        E = E_CAFFE
#        isbad = np.isnan(E)
#        triang.set_mask( np.any(np.where(isbad[triang.triangles], True, False), axis=1) )
        lvlsE = np.arange(0, lvlsE_max+1e-3, 0.25)
        divnorm = colors.TwoSlopeNorm(vmin=np.amin(lvlsE), vcenter=1, vmax=np.amax(lvlsE))
        kwargs_E = dict(levels=lvlsE, norm=divnorm, extend='max', cmap='PuOr_r')
        h = ax.tricontourf(triang, E, **kwargs_E)
        hcb = plt.colorbar(h, ax=ax, ticks=lvlsE[::4], **self.kwargs_cb)
        hcb.set_label(r'$E$')

    def _plot_dlam(self, ax, triang, lami):
        cmap = cmc.lapaz_r
        lvlslam = np.arange(0, 0.8 +1e-3, 0.1)
        kwargs_dlam = dict(levels=lvlslam, extend='max', cmap=cmap)
        h = ax.tricontourf(triang, abs(lami[:,0] - lami[:,1]), **kwargs_dlam)
        hcb = plt.colorbar(h, ax=ax, ticks=lvlslam[::2], **self.kwargs_cb)
        hcb.set_label(r'$\Delta\lambda$')


    def plot_biases(self, FP):
    
        ### Bias analysis 
        
        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        coords,cells, ux,uy,umag,epsE, S,B,H,mask = self.get_np_inputs()
        fabric = self.get_FE_solution(FP, mesh=mesh, boundaries=boundaries)
        u, *_  = self.get_FE_inputs(Q=Q, Q2=Q2, V=V)
        triang = self.get_triang(coords, cells, scale=self.lenscale)

        shearfrac_df, *_ = fabric.enhancementfactor.shearfrac_SSA(u)
        shearfrac        = shearfrac_df.compute_vertex_values(mesh)
#        print(shearfrac[shearfrac>1]) # debug: should be empty list!
        
        chi_df = fabric.enhancementfactor.chi(u, fabric.s)
        chi    = chi_df.compute_vertex_values(mesh)
#        print(chi[chi>1]) # debug: should be empty list!

        ax1,ax2 = self.makefig('', '')

        ax = ax1
#        isbad = np.isnan(shearfrac)
#        triang.set_mask( np.any(np.where(isbad[triang.triangles], True, False), axis=1) )
        h = ax.tricontourf(triang, shearfrac, levels=np.arange(0,1+1e-3,0.1), extend='neither', cmap='Spectral_r')
        hcb = plt.colorbar(h, ax=ax, **self.kwargs_cb)
        #hcb.set_label(r'$\gamma$')
        hcb.set_label(r'{\fontsize{10}{10}\selectfont $\leftarrow$ stretching} \;\; $\gamma$\;\; {\fontsize{10}{10}\selectfont shearing $\rightarrow$}')
        self.plot_domain(ax, triang, mask, scale=self.lenscale, c_floating='k', legend=True)
        
        ax = ax2
#        isbad = np.isnan(chi)
#        triang.set_mask( np.any(np.where(isbad[triang.triangles], True, False), axis=1) )
        h = ax.tricontourf(triang, chi, levels=np.arange(0,1+1e-3,0.1), extend='neither', cmap='PiYG')
        hcb = plt.colorbar(h, ax=ax, **self.kwargs_cb)
        #hcb.set_label(r'$\chi$')
        hcb.set_label(r'{\fontsize{10}{10}\selectfont $\leftarrow$ incompatible} \;\; $\chi$\;\; {\fontsize{10}{10}\selectfont compatible $\rightarrow$}')
        self.plot_domain(ax, triang, mask, scale=self.lenscale, c_floating='k')

        self.savefig(3 if FP=='LROT' else 4)
        
    def plot_domain(self, ax, triang, mask, scale=1, c_floating='limegreen', mesh=False, legend=False, **kwargs):
        if mesh == True: self.plot_mesh(ax, triang)
        self.plot_boundaries(ax, scale=scale, **kwargs)
        self.plot_floating(ax, triang, mask, scale=scale, c=c_floating)
        if legend: 
            ax.add_artist(ax.legend(loc=1, ncols=3, fancybox=False, frameon=False, **self.kwargs_dom))
        
#    def plot_stream(self, ):
#        velcutoff = velscale*10 # don't draw streamlines when velocity is lower than this threshold (m/yr)
#        umag = np.sqrt(np.power(UX,2)+np.power(UY,2))
#        UX[umag < velcutoff] = np.nan # no stream lines over low-velocity areas to avoid noisy plot
#        kwargs_streamplot = dict(density=1.3, linewidth=0.8, color='0.25', broken_streamlines=True)
#        ax.streamplot(X, Y, UX, UY, **kwargs_streamplot)
            
    def plot_CPOs(self, ax1, ax2, fab, lvlmax=0.4, ROTATE_TO_XY=True,
            axsize=0.1, xy=[(0,0),], axloc=[(0,0),], lbl_bbox=(0,0), cnum='w', 
        ):
        
        geo, prj = sfplt.getprojection(rotation=-90, inclination=(0 if ROTATE_TO_XY else 90))
        
        for ii, (x,y) in enumerate(xy):
            W = H = axsize # ax width
            axin = plt.axes([axloc[ii][0]-W/2, axloc[ii][1]-H/2, W, H], projection=prj) 
            axin.set_global()
            nlm = fab.get_nlm(x/self.lenscale, y/self.lenscale)/np.sqrt(4*np.pi) 
            if ROTATE_TO_XY: nlm = sf.rotate_nlm_xz2xy(nlm) # rotate to usual x--y map plane view for SSA models
            sfplt.plotODF(nlm, fab.lm, axin, lvlset=(np.linspace(0.0, lvlmax, 8), lambda x,p:'%.1f'%x), showcb=False, nchunk=None)
            sfplt.plotcoordaxes(axin, geo, axislabels='vuxi', negaxes=False, color='k', fontsize=FS+1.5)
            
            # Markers on map
            num = "%i"%(1+ii)
            peff = [pe.withStroke(linewidth=1.5, foreground='k')]
            axin.add_artist(AnchoredText(num, loc=2, prop=dict(color='k', size=FS+1.5), frameon=False, pad=0.3, bbox_to_anchor=lbl_bbox, bbox_transform=axin.transAxes))
            for jj, ax in enumerate((ax1,ax2)): 
                ax.text(x,y, num, fontsize=FS+3, color='w', path_effects=peff, ha='center', va='center') 

