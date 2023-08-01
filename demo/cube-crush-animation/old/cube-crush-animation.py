#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

import copy, sys, os, code # code.interact(local=locals())
import numpy as np

from scipy import interpolate
import scipy.special as sp

sys.path.insert(0, '../..')
from header import *
from specfabpy import specfabpy as sf 
from sfconstants import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, rc, patches

import warnings
warnings.filterwarnings("ignore")

# plot colors (colorbrewer)
colorax = c_dred
colorbg = "#f7f6ee"

### Fabric time integration 

Nc = 120 # Total number of integration steps taken
Nt0, Nt1 = Nc, int(Nc*0.8)
Nt = Nt0 + Nt1
Nts = 15 # still frames

### Generate fabric time series

class SyntheticFabric():
    
    def __init__(self, L=16): 

        self.L = L # spectral truncation 
        self.lm, self.nlm_len = sf.init(self.L) 
        self.nlm_dummy = np.zeros((self.nlm_len), dtype=np.complex64) 
            
    def M_REG(self, D):

        if self.L <= 8:
            M_REG = sf.M_REG(self.nlm_dummy, D)        
        else:
            nu0, expo = 10, 2
            M_REG = M_REG_custom(nu0, expo, D, sf) # from header.py
    
        return M_REG 
    
    def make_frames(self, dt=25, stressdir='x'):
        
        if stressdir == 'x': stressax, stressaxperp = 0, 2
        if stressdir == 'z': stressax, stressaxperp = 2, 0

        ### Determine dt from Nc and parcel height target
        
        epsii_target = -0.90 # Target: 10% of initial parcel height
        r = 0 # uniaxial compression
        T = 1 # e-folding time scale
        dt = sf.pureshear_strainii_to_t(epsii_target, T)/Nc # time step size for thresshold strain epsii_target in "Nc" time steps 

        ### Construct strain-rate and spin tensor histories
        
        D = np.zeros((Nt+1,3,3)) # strain-rate
        W = np.zeros((Nt+1,3,3)) # spin
        
        # Parcel deformation for plotting
        xyz0_init = [1,1,1]
        xyz0      = np.ones((Nt+1,3))
        xyz0_alt  = np.ones((Nt+1,3))
        strainvec = np.zeros(Nt+1)
        
        # Determine strain history, strain-rate and spin tensors (D and W)
        for ii in np.arange(Nt0+1):
            t = ii*(dt*1)
            ii_ = ii # Nt0-ii
            D[ii_,:,:], W[ii_,:,:] = sf.ugrad_to_D_and_W(sf.pureshear_ugrad(stressax,r,+T))
            Ft = sf.pureshear_F(stressax, r, +T, t)
            strainvec[ii_] = sf.F_to_strain(Ft)[stressax,stressax] # cumulative strain
            xyz0[ii_,:] = np.matmul(Ft, xyz0_init) # parcel side lengths

        for ii in np.arange(Nt1+1):
            t = ii*(dt*0.58)
            ii_ = ii+Nt0
            D[ii_,:,:], W[ii_,:,:] = sf.ugrad_to_D_and_W(sf.pureshear_ugrad(stressax,r,-T))
            Ft = sf.pureshear_F(stressax, r, -T, t)
            strainvec[ii_] = sf.F_to_strain(Ft)[stressax,stressax] # cumulative strain
            xyz0[ii_,:] = np.matmul(Ft, xyz0_init) # parcel side lengths

        ### Fabric evolution 
        
        # Both compression and extension fabric time-series are saved in the same array (nlm_list)
        # Nt0 is the last entry for the first simulation
        
        nlm_list = np.zeros((Nt+1,self.nlm_len), dtype=np.complex64) # The expansion coefficients
        nlm_list[:,0] = 1/np.sqrt(4*np.pi) # Normalized such that N(t=0) = 1 (normalization doesn't matter since we are plotting ODFs)
        eigvals = np.zeros((Nt+1,3))
        
        Eij = np.ones((Nt+1,6))
        e1,e2,e3 = [1,0,0], [0,1,0], [0,0,1]
        (Eij_grain, alpha, n_grain) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
        
        nlm_prev = nlm_list[Nt,:].copy()
        for ii in np.arange(Nt0+1):
            tt = ii
            D_, W_ = D[tt,:,:], W[tt,:,:]
            M = sf.M_LROT(nlm_prev, D_, W_, 1, 0) # Lattice rotation
            M += self.M_REG(D_) # Regularization
            nlm_list[tt,:] = nlm_prev + dt * np.matmul(M, nlm_prev)
            nlm_prev = nlm_list[tt,:]
            Eij[tt,:] = sf.Eij_tranisotropic(nlm_prev, e1,e2,e3, Eij_grain,alpha,n_grain) 
            _,_,_, eigvals[tt,:] = sf.frame(nlm_prev, 'e')
        
        nlm_prev = nlm_list[Nt,:].copy()
        for ii in np.arange(Nt1+1):
            tt = Nt0+ii
            D_, W_ = D[tt,:,:], W[tt,:,:]
            M = sf.M_LROT(nlm_prev, D_, W_, 1, 0) # Lattice rotation
            nu0 = 4
            M += nu0 * sf.M_REG(nlm_prev, D_) # Regularization
            nlm_list[tt,:] = nlm_prev + dt * np.matmul(M, nlm_prev)
            nlm_prev = nlm_list[tt,:]
            Eij[tt,:] = sf.Eij_tranisotropic(nlm_prev, e1,e2,e3, Eij_grain,alpha,n_grain) 
            _,_,_, eigvals[tt,:] = sf.frame(nlm_prev, 'e')
        
        ### Plot

        if 1:
            m = 1
            steps = np.hstack((np.arange(0,Nt0,m),np.arange(Nt0,Nt+1,m))) 
        else: # debug
            steps = [0, int(Nt/2), Nt]
            steps = [0, Nt0, Nt]
        
        for ii in steps:
            
#            if ii <= Nt0: continue # skip comperessive
#            if ii != steps[-1]: continue # skip to last step
            
            scale = 0.4
            fig = plt.figure(figsize=(13*scale,10*scale), facecolor=colorbg)
            
            gs_master = gridspec.GridSpec(1, 2, width_ratios=[0.8,1])
            gs_master.update(top=0.91, bottom=0.13, left=0.095, right=1-0.05, hspace=-0.05, wspace=-0.075)
            
            gsleft = gs_master[0,0].subgridspec(2, 1,  height_ratios=[1,0.3]) 
            axE = fig.add_subplot(gsleft[0, 0])
            
            gs = gs_master[0,1].subgridspec(2, 1,  height_ratios=[1,0.8], hspace=-0.05)
            axParcel = fig.add_subplot(gs[0, 0], projection='3d')
            
            inclination = 55 # view angle
            rot0 = 1 * -90
            rot = 1*-35 + rot0 #*1.4 # view angle
            prj = ccrs.Orthographic(rot, 90-inclination)
            geo = ccrs.Geodetic()            
            axODF = fig.add_subplot(gs[1, 0], projection=prj)
            axODF.set_global()
            
            ### Plot enhancement factors

            xlims = [-1,2]
            ylims = [0,2.5]            
            rect1 = plt.Rectangle((xlims[0],1), np.diff(xlims)[0], ylims[1]-1, color='#deebf7')
            rect2 = plt.Rectangle((xlims[0],0), np.diff(xlims)[0], 1, color='#fee0d2')
            axE.add_patch(rect1)
            axE.add_patch(rect2)
 
            lw = 1.5
            axE.plot(strainvec[[ii,ii]], ylims, '-', c='0.5', lw=lw)
            I = np.concatenate((np.arange(Nt0-1,-1,-1), np.arange(Nt0,Nt+1,1)))

            if stressdir == 'x':
                # Eij = (Exx,Eyy,Ezz,Eyz,Exz,Exy)
                axE.plot(strainvec[I], Eij[I,0], '--k', lw=lw, label=r'$E_{xx}$')
                axE.plot(strainvec[I], Eij[I,2],  ':k', lw=lw, label=r'$E_{zz}$')
                axE.plot(strainvec[I], Eij[I,4], '-k',  lw=lw, label=r'$E_{xz}$')
                axE.plot(strainvec[I], Eij[I,3], '-.k', lw=lw, label=r'$E_{yz}$')
            else:
                raise ValueError('stressdir != x not supported')
                
            axE.text(0.02, 1.50, '{\\bf Softer than}\n{\\bf isotropy}', fontsize=FS-2.5, color='#08519c', ha='center', ma='center')
            axE.text(0.02, 0.20, '{\\bf Harder than}\n{\\bf isotropy}', fontsize=FS-2.5, color='#a50f15', ha='center', ma='center')

            dx=1
            axE.set_xticks(np.arange(xlims[0],xlims[1]+dx, dx))
            axE.set_xticks(np.arange(xlims[0],xlims[1]+dx, dx/2), minor=True)
            axE.set_xlim(xlims)

            dy=1
            axE.set_yticks(np.arange(ylims[0],ylims[1]+dy,dy))
            axE.set_yticks(np.arange(ylims[0],ylims[1]+dy,dy/4), minor=True)
            axE.set_ylim(ylims)
            
            axE.set_xlabel(r'$\epsilon_{zz}$' if stressdir == 'z' else r'$\epsilon_{xx}$')
            axE.set_ylabel('$E_{ij}$')
            axE.set_title('Enhancement factors', pad=8, fontsize=FS-0.5)
            #
            legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'handlelength':1.9, 'labelspacing':0.2}
            hleg = axE.legend(loc=1, fontsize=FS-0.5, **legkwargs)
            hleg.get_frame().set_linewidth(0.7);
            
            ### Plot annotations
            
            x0, y0, dy = -0.125, -0.38, 0.09
            axE.text(x0, y0, r'{\bf specfab}', transform=axE.transAxes, fontsize=FS, horizontalalignment='left')
            axE.text(x0, y0-1*dy, r'github.com/nicholasmr/specfab', transform=axE.transAxes, fontsize=FS-1.5, horizontalalignment='left')
            axE.text(x0, y0-2.0*dy, r'Lattice rotation demo', transform=axE.transAxes, fontsize=FS-1.5, horizontalalignment='left')
            
            x0, y0, dy = +1.875, -0.15, 0.08
            axE.text(x0,y0+2*dy, r'$\lambda_1 = %.2f$'%(eigvals[ii,0]), transform=axE.transAxes, fontsize=FS-1, ha='left')
            axE.text(x0,y0+1*dy, r'$\lambda_2 = %.2f$'%(eigvals[ii,1]), transform=axE.transAxes, fontsize=FS-1, ha='left')
            axE.text(x0,y0+0*dy, r'$\lambda_3 = %.2f$'%(eigvals[ii,2]), transform=axE.transAxes, fontsize=FS-1, ha='left')
            
            ### Plot parcel geometry
            
            plot_parcel(axParcel, xyz0[ii,:], 0,0,0, color='0.5', colorax=colorax, scale=0.60, plotaxlbls=True) # see header.py

            ### Plot ODF

            plot_ODF(nlm_list[ii,:], self.lm, ax=axODF, rot0=rot0, cblabel='ODF') # see header.py
            
            axParcel.set_facecolor(colorbg)
            axODF.text(rot0-85, 82, r'$\vu{z}$', c=colorax, horizontalalignment='left', transform=geo)
            axODF.text(rot0-90-8*1, -8, r'$\vu{x}$', c=colorax, horizontalalignment='left', transform=geo)
            axODF.text(rot0-5*1, -5, r'$\vu{y}$', c=colorax, horizontalalignment='left', transform=geo)
            
            ### Save frame in correct order (file name number) for generating animation below

            if ii <= Nt0:
                fout = 'frames/frame%04i.png'%(ii+Nts)
                plt.savefig(fout, dpi=300)

            if ii > Nt0:
                fout = 'frames/frame%04i.png'%(ii+Nts*3-1)
                plt.savefig(fout, dpi=300)

            print('Saving %i :: %s'%(ii,fout))    
            plt.close()

    
### Generate frames and make animation

os.system('mkdir -p frames')    
synfab = SyntheticFabric()
synfab.make_frames()    

# Still frames
k = Nt0+2*Nts
for ii in np.concatenate((np.arange(0,Nts), k+np.arange(0,Nts))):
    os.system('cp frames/frame%04d.png frames/frame%04d.png'%(Nts, ii))

for ii in np.arange(0,Nts+1): os.system('cp frames/frame%04d.png frames/frame%04d.png'%(Nt0+Nts-1, Nt0+Nts+ii))
for ii in np.arange(0,Nts+1): os.system('cp frames/frame%04d.png frames/frame%04d.png'%(Nt+3*Nts-1, Nt+3*Nts+ii))

# Make .avi and GIF?
if 1:
    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/frame%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p cube-crush.avi')

    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/frame%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p cube-crush-for-gif.avi')
    os.system('rm cube-crush.gif')
    os.system('ffmpeg -i cube-crush-for-gif.avi -vf "fps=22,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 cube-crush.gif')
    os.system('rm cube-crush-for-gif.avi')    
       
