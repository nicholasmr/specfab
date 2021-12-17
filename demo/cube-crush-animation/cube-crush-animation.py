#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021

import copy, sys, os, code # code.interact(local=locals())
import numpy as np

from scipy import interpolate
import scipy.special as sp

sys.path.insert(0, '..')
from header import *
from specfabpy import specfabpy as sf # requires the spectral fabric module to be compiled!

s2yr   = 3.16887646e-8
yr2s   = 31556926    
yr2kyr = 1e-3 

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, rc, patches

# plot colors (colorbrewer)
colorax = c_dred

class SyntheticFabric():
    
    def __init__(self, L=16, nu=3e-3): 

        self.L  = L # spectral truncation 
        self.nu0 = nu # Regularization strength (calibrated for L and strain-rate)
               
        self.nlm_len = sf.init(self.L) 
        self.lm = sf.get_lm(self.nlm_len)
    
    ####################################

    def make_profile(self, dt=25):
        
        STRESSDIRECTION = 'x'
        #STRESSDIRECTION = 'z'
        
        if STRESSDIRECTION == 'x': 
            STRESSAX = 0
            STRESSAXPERP = 2
            
        if STRESSDIRECTION == 'z': 
            STRESSAX = 2
            STRESSAXPERP = 0
    
        deformExpr1 = {'ax':STRESSDIRECTION, 't_c':+1000, 'r':0, 'Gamma0':0}
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
            strain_zz[Nt0-ii] = PS.strain(t)[STRESSAX,STRESSAX]
            xyz0[Nt0-ii,:] = np.matmul(PS.F(t), xyz0_init)
            W[Nt0-ii,:,:], D[Nt0-ii,:,:] = PS.W(), PS.D()            

        for ii in np.arange(Nt1+1):
            t = ii*dt
            PS = PureShear(-deformExpr['t_c']*yr2s, deformExpr['r'], ax=deformExpr['ax'])
            strain_zz[Nt0+ii] = PS.strain(t)[STRESSAX,STRESSAX]
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
            inclination = 55 # view angle
            rot0 = 1 * -90
            rot = 1*-35 + rot0 #*1.4 # view angle
            prj = ccrs.Orthographic(rot, 90-inclination)
            geo = ccrs.Geodetic()            
            axODF = fig.add_subplot(gs[1, 1], projection=prj)
            axODF.set_global()
            
            ##
 
            lw = 1.5
            axE.plot(Eij[:,STRESSAX,STRESSAX], strain_zz, '--k',  lw=lw, label=r'$E_{zz}$' if STRESSDIRECTION == 'z' else r'$E_{xx}$')
            axE.plot(Eij[:,STRESSAXPERP,STRESSAX], strain_zz, '-k', lw=lw, label=r'$E_{xz}$' if STRESSDIRECTION == 'z' else r'$E_{zx}$')

            ylims = [-1,2]
            dy=0.5
            axE.set_yticks(np.arange(ylims[0],ylims[1]+dy, dy))
            axE.set_yticks(np.arange(ylims[0],ylims[1]+dy, dy/2), minor=True)
            axE.set_ylim(ylims)

            xlims = [0,2.5]
            axE.plot(xlims,strain_zz[[ii,ii]], ':', c=c_green, lw=lw)
            dx=1
            axE.set_xticks(np.arange(xlims[0],xlims[1]+dx,dx))
            axE.set_xticks(np.arange(xlims[0],xlims[1]+dx,dx/4), minor=True)
            axE.set_xlim(xlims)
            
            axE.set_ylabel(r'$\epsilon_{zz}$' if STRESSDIRECTION == 'z' else r'$\epsilon_{xx}$')
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
            
            plot_parcel(axParcel, xyz0[ii,:], 0,0,0, color='0.5', colorax=colorax, scale=0.60, plotaxlbls=True)
            plot_ODF(nlm_list[ii,:], self.lm, ax=axODF, rot0=rot0, cblabel='ODF')        
        
            mrk='.'
            ms=5
            axODF.plot([0],[90],mrk, ms=ms, c=colorax, transform=geo)
            axODF.plot([rot0-90],[0],mrk, ms=ms,  c=colorax, transform=geo)
            axODF.plot([rot0],[0],mrk, ms=ms,  c=colorax, transform=geo)
            axODF.text(rot0-80, 70, r'$\vu{z}$', c=colorax, horizontalalignment='left', transform=geo)
            axODF.text(rot0-90+8, -5, r'$\vu{x}$', c=colorax, horizontalalignment='left', transform=geo)
            axODF.text(rot0-15, -5, r'$\vu{y}$', c=colorax, horizontalalignment='left', transform=geo)
            
            ###

            print('Saving %i'%(ii))        
            if ii <= Nt0:
                fout1 = 'frames/frame%04i.png'%(Nt0-ii)
                fout2 = 'frames/frame%04i.png'%(Nt0+ii)
                plt.savefig(fout1, dpi=300)
                os.system('cp %s %s'%(fout1,fout2))

            if ii > Nt0:
                fout1 = 'frames/frame%04i.png'%(Nt0 +ii)
                fout2 = 'frames/frame%04i.png'%(2*Nt - (ii-Nt0))
                plt.savefig(fout1, dpi=300)
                print('2*Nt=%i :: %i :: %i'%(2*Nt, Nt0 +ii, 2*Nt - (ii-Nt0)))
                os.system('cp %s %s'%(fout1,fout2))

            plt.close()
    
###

os.system('mkdir -p frames')    
synthfab = SyntheticFabric()
synthfab.make_profile()    
os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 1 -i frames/frame%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p cube-crush.avi')

# Make GIF
if 0:
    os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/frame%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p cube-crush-for-gif.avi')
    os.system('ffmpeg -i cube-crush-for-gif.avi -vf "fps=25,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 cube-crush.gif')
    os.system('rm cube-crush-for-gif.avi')
    
       
