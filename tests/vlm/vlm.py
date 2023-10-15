# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

"""
Given two ODFs of orthogonal axes (coefficieints blm and nlm), calculate ODF of the third axis (coefficients vlm)
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

#----------------------
# Idealized structure tensors for test cases
#----------------------

x,y,z = [1,0,0],[0,1,0],[0,0,1]

L = 2
x2 = sf.nlm_ideal(x,0,L) 
y2 = sf.nlm_ideal(y,0,L)
z2 = sf.nlm_ideal(z,0,L)
G2xy = sf.nlm_ideal(z,np.pi/2,L)
G2yz = sf.nlm_ideal(x,np.pi/2,L)

L = 4
x4 = sf.nlm_ideal(x,0,L) 
y4 = sf.nlm_ideal(y,0,L)
z4 = sf.nlm_ideal(z,0,L)
G4xy = sf.nlm_ideal(z,np.pi/2,L)
G4yz = sf.nlm_ideal(x,np.pi/2,L)

a2_pairs = ((x2,y2), (y2,z2), (G2xy,z2), (G2yz,x2))
a4_pairs = ((x4,y4), (y4,z4), (G4xy,z4), (G4yz,x4))

#----------------------
# Setup figure
#----------------------

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

dpi, scale = 200, 2.5
figsize=(2*scale,2*1.5*scale)
a = 0.00
kwargs_gsupdate = {'left':a, 'right':1-a, 'top':0.99, 'bottom':0.07, 'wspace':0.015*10, 'hspace':0.6}

geo, prj = sfplt.getprojection(rotation=45, inclination=45)

lvlset = (np.linspace(0,0.4,9), lambda x,p:'%.1f'%x) # custom level set: (list of levels, how to format colorbar tick labels)

def make_fig(rows=4,cols=3,figsize=figsize):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(rows,cols)
    gs.update(**kwargs_gsupdate)
    return fig, gs

def get_axes(gs,ii):
    ax1 = plt.subplot(gs[ii, 0], projection=prj)
    ax2 = plt.subplot(gs[ii, 1], projection=prj)
    ax3 = plt.subplot(gs[ii, 2], projection=prj)
    ax1.set_global(); ax2.set_global(); ax3.set_global()
    return ax1,ax2,ax3

#----------------------
# Test L=2 truncation
#----------------------

if 1:

    print('Calculating for L=2')
    lm, nlm_len = sf.init(2)
    fig, gs = make_fig()
        
    for ii, (blm, nlm) in enumerate(a2_pairs):

        (ax1,ax2,ax3) = get_axes(gs,ii)

        vlm = sf.a2_to_nlm(sf.a2_orth(blm,nlm))

        if 0:        
            print('------ a2 for blm, nlm, vlm ---------')        
            print(sf.a2(blm))
            print(sf.a2(nlm))
            print(sf.a2_orth(blm,nlm))
            print('---------------')

        sfplt.plotODF(blm, lm, ax1, cblabel=r'$b/N$ ($L=2$)', lvlset=lvlset)
        sfplt.plotcoordaxes(ax1, geo, axislabels='vuxi')

        sfplt.plotODF(nlm, lm, ax2, cblabel=r'$n/N$ ($L=2$)', lvlset=lvlset)
        sfplt.plotcoordaxes(ax2, geo, axislabels='vuxi')
        
        sfplt.plotODF(vlm, lm, ax3, cblabel=r'$v/N$ ($L=2$)', lvlset=lvlset)
        sfplt.plotcoordaxes(ax3, geo, axislabels='vuxi')

    fout = 'vlm-L2-deltafunc.png'
    print('Saving %s\n'%(fout))
    plt.savefig(fout, dpi=dpi)

#----------------------
# Test L=4 truncation
#----------------------

if 1:

    print('Calculating for L=4')
    lm, nlm_len = sf.init(4)
    fig, gs = make_fig()
        
    for ii, (blm, nlm) in enumerate(a4_pairs):

        (ax1,ax2,ax3) = get_axes(gs,ii)
        
        vlm = sf.a4_to_nlm(sf.a4_orth(blm,nlm))

        sfplt.plotODF(blm, lm, ax1, cblabel=r'$b/N$ ($L=4$)', lvlset=lvlset)
        sfplt.plotcoordaxes(ax1, geo, axislabels='vuxi')

        sfplt.plotODF(nlm, lm, ax2, cblabel=r'$n/N$ ($L=4$)', lvlset=lvlset)
        sfplt.plotcoordaxes(ax2, geo, axislabels='vuxi')
        
        sfplt.plotODF(vlm, lm, ax3, cblabel=r'$v/N$ ($L=4$)', lvlset=lvlset)
        sfplt.plotcoordaxes(ax3, geo, axislabels='vuxi')

    fout = 'vlm-L4-deltafunc.png'
    print('Saving %s\n'%(fout))
    plt.savefig(fout, dpi=dpi)
    
#----------------------
# Test with specfab under shear
#----------------------

if 1:

    for Ltrunc in [2,4]:

        lm_L, len_L = sf.init(Ltrunc) # just for saving returns

        L=12
        print('Calculating for L=%i, Ltrunc=%i for evolving fabric'%(L,Ltrunc))

        exprtypes = ['ss','uc']
        
        for exprtype in exprtypes:

            ### Velocity gradient tensor experienced by parcel
            if exprtype=='ss':
                pl, T = 1, 1
                D, W = sf.ugrad_to_D_and_W(sf.simpleshear_ugrad(pl, T))
                Nt = 6*25 # Number of time steps
            if exprtype=='uc':
                ax, r, T = 2, 0, 1
                D, W = sf.ugrad_to_D_and_W(sf.pureshear_ugrad(ax, r, T))
                Nt = 4*25 # Number of time steps

            ### Numerics 
            dt = 1*0.025 # Time-step size

            ### Initialize fabric as isotropic
            lm, nlm_len = sf.init(L) 
            nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # Array of expansion coefficients
            nlm[0,0] = 1/np.sqrt(4*np.pi) # Normalized ODF at t=0
            blm = nlm.copy()
            vlm = nlm.copy()
            nlm_iso = nlm[0,:].copy()

            #### Euler integration of lattice rotation + regularization
            for tt in np.arange(1,Nt):
                
                blm_prev    = blm[tt-1,:] # Previous solution
                nlm_prev    = nlm[tt-1,:] # Previous solution
                
                M_LROT_b = sf.M_LROT(nlm_prev, D, W, -1, 0)
                M_LROT_n = sf.M_LROT(nlm_prev, D, W, +1, 0)
                M_REG  = sf.M_REG(nlm_prev, D)
                
                blm[tt,:] = blm_prev + dt*np.matmul(M_LROT_b+M_REG, blm_prev) # Euler step
                nlm[tt,:] = nlm_prev + dt*np.matmul(M_LROT_n+M_REG, nlm_prev) # Euler step

                if Ltrunc == 2: vlm[tt,:sf.L2len] = sf.a2_to_nlm(sf.a2_orth(blm[tt,:],nlm[tt,:]))
                if Ltrunc == 4: vlm[tt,:sf.L4len] = sf.a4_to_nlm(sf.a4_orth(blm[tt,:],nlm[tt,:]))

            fig, gs = make_fig(rows=2, figsize=(figsize[0],0.4*figsize[1]))
            (ax1,ax2,ax3) = get_axes(gs,0)
                
            tt = -1

            sfplt.plotODF(blm[tt,:len_L], lm, ax1, cblabel=r'$b/N$ ($L=%i$)'%(Ltrunc), lvlset=lvlset)
            sfplt.plotcoordaxes(ax1, geo, axislabels='vuxi')

            sfplt.plotODF(nlm[tt,:len_L], lm, ax2, cblabel=r'$n/N$ ($L=%i$)'%(Ltrunc), lvlset=lvlset)
            sfplt.plotcoordaxes(ax2, geo, axislabels='vuxi')
            
            sfplt.plotODF(vlm[tt,:len_L], lm, ax3, cblabel=r'$v/N$ ($L=%i$)'%(Ltrunc), lvlset=lvlset)
            sfplt.plotcoordaxes(ax3, geo, axislabels='vuxi')

            if 1:
            
                print('\n...and test for two isotropic distributions')

                if Ltrunc==2: 
                    a2_v  = sf.a2_orth(nlm_iso,nlm_iso)
                    vlm = sf.a2_to_nlm(a2_v)
                    print('a2_v for b,n=iso:\n',a2_v)
                
                if Ltrunc==4: 
                    a4_v  = sf.a4_orth(nlm_iso,nlm_iso)
                    vlm = sf.a4_to_nlm(a4_v)
        #            print('a4_v for b,n=iso:\n',a4_v)

                print('vlm (should be isotropic) = ', vlm)

                (ax1,ax2,ax3) = get_axes(gs,1)
                sfplt.plotODF(vlm[:len_L], lm_L, ax3, cblabel=r'$v/N$ | isotropic CPO', lvlset=lvlset)
                sfplt.plotcoordaxes(ax3, geo, axislabels='vuxi')

            fout = 'vlm-L%i-sf-%s.png'%(Ltrunc,exprtype)
            print('Saving %s'%(fout))
            plt.savefig(fout, dpi=dpi)
            
