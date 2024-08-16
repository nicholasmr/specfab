# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

"""
Script for processing experimental data that generates a2 and a4 from measured grain orientations
"""

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import quaternion as qt # pip3 install numpy-quatern
import pandas as pd
import pickle, glob

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

sys.path.insert(0, '..')
from localheader import * 
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)

PLOT_OVERVIEW_FIGURE = True

#--------------------
# Select experiments to process
#--------------------

# See definitions in experiments.py 

# All data
#experiments = (expr_Kim2022, expr_Yabe2020, expr_Miyazaki2013_uc, expr_Miyazaki2013_ue, expr_Bernard2019, expr_Kumamoto2019)
experiments = (expr_Yabe2020, expr_Miyazaki2013_uc, expr_Miyazaki2013_ue, expr_Bernard2019, expr_Kumamoto2019)

# Subsets
#experiments = (expr_Yabe2020, expr_Miyazaki2013_uc, expr_Miyazaki2013_ue, )
#experiments = (expr_Bernard2019, expr_Kumamoto2019,expr_Kim2022)
#experiments = (expr_Miyazaki2013_uc,expr_Miyazaki2013_ue)
experiments = (expr_Boneh2014, )

#experiments = (expr_DREX_ss,)

#--------------------
# Process selected experiments
#--------------------

for expr in experiments:

    corr = {'b20':[],'b40':[],'b60':[], 'n20':[],'n40':[],'n60':[], 'F':np.zeros((3,3))} # container for correlation data
    
    print('\n----------------------------------')
    print('Processing experiments in %s/'%(expr['path']))
    print('----------------------------------')

    flist = expr['flist']
    expr_fnames = [os.path.basename(f)[0:] for f in flist] 

    for ff, fname0 in enumerate(flist):

        if expr['source'] == 'drex':
            fname = '%s-F.csv'%(fname0[:-4])
            fullfname = 'data/%s/%s'%(expr['path'],fname)
            df = pd.read_csv(fullfname, header=None, sep=',') 
            corr['F'] = df.to_numpy()

        for ii in [1,2,3]:
        
            if ii==1: axname = 'b' 
            if ii==2: axname = 'n'
            if ii==3: axname = 'v'            

            ### Load sample
            
            # -1 => smax, 0 => girdle
            Ilami = None
            if axname == 'b':
                if expr['type'] == 'ss': Ilam1 = -1 
                if expr['type'] == 'cc': Ilam1 = -1
                if expr['type'] == 'uc': Ilam1 = 0
                if expr['type'] == 'ue': Ilam1 = -1
            if axname == 'n':
                if expr['type'] == 'ss': Ilam1 = -1
                if expr['type'] == 'cc': Ilam1 = -1
                if expr['type'] == 'uc': Ilam1 = -1
                if expr['type'] == 'ue': Ilam1 = 0

            fname = '%s-%i.csv'%(fname0, ii)
            (q, qr, m, mnew, caxes, nlm, nlmr, lm) = load_sample(fname, expr, Ilam1=Ilam1)
            (qlat, qcolat, qlon) = q
            (qlatr, qcolatr, qlonr) = qr
            (m_lat, m_colat, m_lon) = m
            (mnew_lat, mnew_colat, mnew_lon) = mnew

            ### Lowest-order representation for reference

            nlm_L2 = nlm.copy()
            nlm_L2[sf.L2len:] = 0
            
            nlmr_L2 = nlmr.copy()
            nlmr_L2[sf.L2len:] = 0
            
            ### Save correlation for later plotting
            
            if axname == 'b':
                corr['b20'].append(np.real(nlmr[sf.I20])) # nl0 are identically real 
                corr['b40'].append(np.real(nlmr[sf.I40]))
                corr['b60'].append(np.real(nlmr[sf.I60]))
            
            if axname == 'n':
                corr['n20'].append(np.real(nlmr[sf.I20])) # nl0 are identically real 
                corr['n40'].append(np.real(nlmr[sf.I40]))
                corr['n60'].append(np.real(nlmr[sf.I60]))
            
            """
            Setup figure for plotting the derived ODF and rotated mi' of this sample.
            This is an overview figure meant for verifying the sample was processed correctly.
            """

            if PLOT_OVERVIEW_FIGURE:
            
                scale = 2.2
                fig = plt.figure(figsize=(2.7*scale,3*scale))
                gs = gridspec.GridSpec(2,3)
                gs.update(left=0.03, right=1-0.03/3, top=0.99, wspace=0.015*18, hspace=0.35)

                geo, prj = sfplt.getprojection(rotation=0, inclination=40)
                ax1,ax2,ax3    = [plt.subplot(gs[0, ii], projection=prj) for ii in range(3)]
                ax1r,ax2r,ax3r = [plt.subplot(gs[1, ii], projection=prj) for ii in range(3)]
                for axi in (ax1,ax2,ax3, ax1r,ax2r,ax3r): axi.set_global()

                kwargs_ODF = dict(lvlset=(np.linspace(0,0.4,9),lambda x,p:'%.1f'%x), cbtickintvl=4, cblabel=r'$%s/N$'%(axname))
                kwargs_dp  = dict(ms=0.5, marker='o', c=sfplt.c_blue, transform=geo) # data points
                kwargs_m   = dict(marker='o', markersize=4, c=sfplt.c_green, zorder=10, transform=geo)
                kwargs_cax = dict(axislabels='vuxi', color=sfplt.c_dred)

                ### Plot -- true frame
                
                sfplt.plotODF(nlm,    lm, ax1, **kwargs_ODF)
                sfplt.plotODF(nlm_L2, lm, ax2, **kwargs_ODF)
                sfplt.plotODF(nlm,    lm, ax3, **kwargs_ODF)

                plot_points(ax1, *get_deg(qlat, qlon), **kwargs_dp)

                for axi in (ax1,ax2,ax3): 
                    sfplt.plotcoordaxes(axi, geo, **kwargs_cax)
                    plot_points(axi, *get_deg(sfdsc.colat2lat(m_colat), m_lon), **kwargs_m)

                ### Plot -- rotated frame
                
                sfplt.plotODF(nlmr,    lm, ax1r, **kwargs_ODF)
                sfplt.plotODF(nlmr_L2, lm, ax2r, **kwargs_ODF)
                sfplt.plotODF(nlmr,    lm, ax3r, **kwargs_ODF)

                plot_points(ax1r, *get_deg(qlatr, qlonr), **kwargs_dp)
                            
                for axi in (ax1r,ax2r,ax3r): 
                    sfplt.plotcoordaxes(axi, geo, **kwargs_cax)
                    plot_points(axi, *get_deg(sfdsc.colat2lat(mnew_colat), mnew_lon), **kwargs_m)
                
                ### Save figure

                fout = r'data/%s/overview-%s.png'%(expr['path'], fname)
                print('- Saving %s'%(fout))
                plt.savefig(fout, dpi=250)
                plt.close()

    ### Save all correlations from the files of this experiment
    fcorr = "observed-states/%s.p"%(expr['path']) # save this experiment's correlations in this file
    print('\n=== Saving correlations: %s ==='%(fcorr))
    pickle.dump(corr, open(fcorr, "wb"))
       
