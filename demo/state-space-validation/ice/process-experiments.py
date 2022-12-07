# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Script for processing experimental data that generates a2 and a4 from measured c-axis ensembles
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import quaternion as qt # pip3 install numpy-quatern
import pandas as pd
import pickle, glob
from progress.bar import Bar

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

from localheader import * 
from experiments import * # experiment definitions (data structures for experiment files, etc.)
sys.path.insert(0, '../../../demo') # for importing local specfabpy build (if available) and common python header
from header import * # contains matplotlib setup etc.
from specfabpy import specfabpy as sf

#--------------------
# Select experiments to process
#--------------------

# See definitions in experiments.py 

# All data
#experiments = (expr_GRIP, expr_LAWDOME,  expr_EGRIP_MAIN, expr_SPICE, expr_EGRIP_S5, expr_Priestly, expr_Qi, expr_Fan_10,expr_Fan_4, expr_Hunter) 

# Subsets
#experiments = ()
#experiments = (expr_Qi, expr_Fan_10,expr_Fan_4, expr_Hunter)
#experiments = (expr_Qi,expr_Fan_30,expr_Fan_20,expr_Fan_10,expr_Fan_4)
#experiments = (expr_Fan_10,expr_Fan_4)
#experiments = (expr_Hunter,)
experiments = (expr_Priestley,)

#--------------------
# Process selected experiments
#--------------------

for expr in experiments:

    corr = {'n20': [], 'n40': []} # container for data correlation 
    
    print('\n----------------------------------')
    print('Processing experiments in %s/'%(expr['path']))
    print('----------------------------------')

    flist = expr['flist']
    expr_fnames = [os.path.basename(f)[0:] for f in flist] # name according to three first letters of file name

    for ff, fname in enumerate(flist):

        fullfname = 'data/%s/%s'%(expr['path'],fname)
        print('\n*** Loading %s'%(fullfname))
        
        df = pd.read_csv(fullfname, skiprows=expr['skiprows'], sep=expr['sep']) # header=expr['headerlen'], 

        # For SPICE files remove rows with "unrecognized grains"
        m = ~df.astype(str).apply(lambda x: x.str.contains('grain', regex=False)).any(axis=1)
        df = df[m]            
        
        if expr['coords'] == 'quaternion': # Processed EBSD data in MTEX (MATLAB) is saved as quaternion coords.
            #qs_comps = df.to_numpy()[:,:-1] # @TODO :: Does not have four components (columns) with this method if grain sizes are also saved in MTEX.
            qs_comps = df.to_numpy()[:,:4] 
            qs = qt.as_quat_array(qs_comps)
            sphcoords = qt.as_spherical_coords(qs)
            qcolat, qlon = sphcoords[:,0], sphcoords[:,1] # COLAT [0;pi], LON [0;2pi]
        #    qcolat, qlon = np.array([np.pi/2 * 1,]), np.array([np.pi/2 * 1,]) # debug
            qlat = np.pi/2 - qcolat 
            
        else: # c-axes given by spherical coordinates?
            sphcoords = df.to_numpy()
            qlat, qlon = np.deg2rad(sphcoords[:,expr['I_lat']].astype(np.float64)), np.deg2rad(sphcoords[:,expr['I_azi']].astype(np.float64))
            qcolat = np.pi/2 - qlat
            if expr['iscolat']:
                qlat, qcolat = qcolat, qlat # swap

        # Add reflected vectors to c-axis bundle; ODFs are antipodally symmetric, so doesn't change statistics, but is useful when plotting a c-axes bundle on top of its ODF.
        if 1: 
            qcolat = np.hstack((qcolat, np.pi-qcolat))
            qlon   = np.hstack((qlon, qlon-np.pi))
            qlat = np.pi/2 - qcolat 
            
        ### Determine ODF from a2 and a4
        
        lm, nlm_len = sf.init(4)
        nlm_L4 = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients
        nlm_L2 = nlm_L4.copy()
        caxes = np.array([ [np.cos(p)*np.sin(t), np.sin(p)*np.sin(t), np.cos(t)] for t, p in zip(qcolat,qlon) ]) # construct c-axis from angles
        a2 = np.array([ np.einsum('i,j',c,c)         for c in caxes]).mean(axis=0) # construct <cc>
        a4 = np.array([ np.einsum('i,j,k,l',c,c,c,c) for c in caxes]).mean(axis=0) # construct <cccc>
        nlm_L2[:6]  = sf.a2_to_nlm(a2) # for comparing with a4 method
        nlm_L4[:16] = sf.a4_to_nlm(a4)

        ### Rotated frame 
           
        Ilami = 0 if expr['type'] == 'ue' else 2 # sym. axis = largest eig dir. for single max, smallest eig dir. for girdle
        (v1_colat, v1_lon, _) = get_v_angles(a2, Ilami=Ilami)
        #v1_colat = 0 # only horiz rotation (for debugging)
        
        print('true v1 colat, lon = %f, %f (deg.) '%(np.rad2deg(v1_colat), np.rad2deg(v1_lon)))
        # THIS is the nlm array from which spectral coefs are derived for correlation
        nlmr_L4 = sf.rotate_nlm(nlm_L4, 0, -v1_lon) # The rotation is the composition of two rotations 
        nlmr_L4 = sf.rotate_nlm(nlmr_L4, +v1_colat, 0) # ... second rotation
        # *_L2 is used only for plotting!
        nlmr_L2 = nlmr_L4.copy()
        nlmr_L2[6:] = 0        
        
        # verify that rotated a2 frame has v1 vertical
        (v1new_colat, v1new_lon, v1new) = get_v_angles(sf.a2(nlmr_L4))
        print('new  v1 colat, lon = %f, %f (deg.) '%(np.rad2deg(v1new_colat), np.rad2deg(v1new_lon)))
        print('new  v1 = ', v1new)
        
        Ry = np.matrix([[np.cos(-v1_colat),0,np.sin(-v1_colat)],[0,1,0],[-np.sin(-v1_colat),0,np.cos(-v1_colat)]]) # R_y
        Rz = np.matrix([[np.cos(-v1_lon),-np.sin(-v1_lon),0],[np.sin(-v1_lon),np.cos(-v1_lon),0],[0,0,1]]) # R_z
        caxes = np.array([np.einsum('ij,jk,k',Ry,Rz,cax) for cax in caxes])
        qlonr   = np.arctan2(caxes[:,1],caxes[:,0]) # arctan2(y,x)
        qcolatr = np.arccos(caxes[:,2])  
        qlatr   =  np.pi/2 - qcolatr
        
        ### Save correlation for later plotting
        
        corr['n20'].append(np.real(nlmr_L4[ 3])) # n_2^0 (is real by def.)
        corr['n40'].append(np.real(nlmr_L4[10])) # n_4^0 (is real by def.)
        
        """
        Setup figure for plotting the derived ODF and rotated c-axes of this sample.
        This is an overview figure meant for verifying the sample was processed correctly.
        """

        if 1:
        
            dpi, scale = 250, 2.2
            fig = plt.figure(figsize=(2.7*scale,3*scale))
            gs = gridspec.GridSpec(2,3)
            a = 0.03
            gs.update(left=a, right=1-a/3, top=0.99, wspace=0.015*18, hspace=0.35)
            gs.update(wspace=0.015*18, hspace=0.35)

            geo = ccrs.Geodetic()
            rot0 = 0
            rot = 1.4*rot0 # view angle
            inclination = 40 # view angle
            prj = ccrs.Orthographic(rot, 90-inclination)
            ax1 = plt.subplot(gs[0, 0], projection=prj)
            ax2 = plt.subplot(gs[0, 1], projection=prj)
            ax3 = plt.subplot(gs[0, 2], projection=prj)
            ax1r = plt.subplot(gs[1, 0], projection=prj)
            ax2r = plt.subplot(gs[1, 1], projection=prj)
            ax3r = plt.subplot(gs[1, 2], projection=prj)

            ax1.set_global(); ax2.set_global(); ax3.set_global() 
            ax1r.set_global(); ax2r.set_global(); ax3r.set_global() 

            ### Plot

            qlatd, qlond = get_deg(qlat, qlon)
            v1latd, v1lond = get_deg(np.pi/2-v1_colat, v1_lon)
            ms = 0.3

            plot_ODF(nlm_L4, lm, ax=ax1, cmap='Greys', cblabel=r'$\psi/N$')
            ax1.plot(qlond, qlatd, ls='none', marker='x', markersize=ms, c='tab:blue', transform=geo) 
            ax1.plot(v1lond, v1latd, ls='none', marker='o', markersize=4, c='tab:green', transform=geo) 
            plot_axes(ax1, geo)
            
            plot_ODF(nlm_L2, lm, ax=ax2, cmap='Greys', cblabel=r'$\psi/N, L=2$')
            ax2.plot(v1lond, v1latd, ls='none', marker='o', markersize=4, c='tab:green', transform=geo) 
            plot_axes(ax2, geo)

            plot_ODF(nlm_L4, lm, ax=ax3, cmap='Greys', cblabel=r'$\psi/N, L=4$')
            ax3.plot(v1lond, v1latd, ls='none', marker='o', markersize=4, c='tab:green', transform=geo) 
            plot_axes(ax3, geo)

            ### Plot again for rotated frame
            
            qlatrd, qlonrd = get_deg(qlatr, qlonr)
            v1newlatd, v1newlond = get_deg(np.pi/2-v1new_colat, v1new_lon)
            
            plot_ODF(nlmr_L4, lm, ax=ax1r, cmap='Greys', cblabel=r'$\psi/N$ (rot. frame)')
            ax1r.plot(qlonrd, qlatrd, ls='none', marker='x', markersize=ms, c='tab:blue', transform=geo) 
            plot_axes(ax1r, geo)
            ax1r.plot(v1newlond, v1newlatd, ls='none', marker='o', markersize=4, c='tab:green', transform=geo) 
            
            plot_ODF(nlmr_L2, lm, ax=ax2r, cmap='Greys', cblabel=r'$\psi/N, L=2$ (rot. frame)')
            plot_axes(ax2r, geo)
            ax2r.plot(v1newlond, v1newlatd, ls='none', marker='o', markersize=4, c='tab:green', transform=geo) 
            
            plot_ODF(nlmr_L4, lm, ax=ax3r, cmap='Greys', cblabel=r'$\psi/N, L=4$ (rot. frame)')
            plot_axes(ax3r, geo)
            ax3r.plot(v1newlond, v1newlatd, ls='none', marker='o', markersize=4, c='tab:green', transform=geo) 
            
            ### Save figure

            fout = r'data/%s/overview-%s.png'%(expr['path'], expr_fnames[ff])
            print('*** Saving %s'%(fout))
            plt.savefig(fout, dpi=dpi)
            plt.close()


    ### Save all correlations from the files of this experiment
    fcorr = "observed-states/%s.p"%(expr['path']) # save this experiment's correlations in this file
    print('\n=== Saving correlations: %s ==='%(fcorr))
    pickle.dump(corr, open(fcorr, "wb"))
       
