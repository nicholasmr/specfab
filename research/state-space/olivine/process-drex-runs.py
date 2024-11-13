# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

"""
Script for processing D-REX experiments
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

#--------------------
# Select experiments to process
#--------------------

# See definitions in experiments.py 

experiments = (expr_DREX_uc, expr_DREX_ue, expr_DREX_ss)
#experiments = (expr_DREX_uc,)

#--------------------
# Process selected experiments
#--------------------

lm, nlm_len = sf.init(6)

USE_DREX_FOR_v = 0

def load_axes(fullfname):
    df = pd.read_csv(fullfname, skiprows=expr['skiprows'], header=None, sep=expr['sep'])
    vj = df.to_numpy()
    return vj
    
def axes_to_nlm(axes):
    a6 = np.array([ np.einsum('i,j,k,l,m,n',r,r,r,r,r,r) for r in axes]).mean(axis=0) # construct <c^6>
    return sf.a6_to_nlm(a6) 

for expr in experiments:

    print('\n----------------------------------')
    print('Processing %s/'%(expr['path']))
    print('----------------------------------')

    flist = expr['flist'] # files of a given deformation experiment

    N = len(flist)
    F = np.zeros((3,3,N))
    nlm_1 = np.zeros((N,nlm_len), dtype=np.complex64)
    nlm_2 = nlm_1.copy()
    nlm_3 = nlm_1.copy()
    
    for ff, fname0 in enumerate(flist):

        print('*** Loading %s splitted files'%(fname0))
        
        # Axis orientations
        fname1 = '%s-%i.csv'%(fname0[:-4],1)
        fname2 = '%s-%i.csv'%(fname0[:-4],2)
        fname3 = '%s-%i.csv'%(fname0[:-4],3)
        axes1 = load_axes('data/%s/%s'%(expr['path'],fname1))
        axes2 = load_axes('data/%s/%s'%(expr['path'],fname2))
        axes3 = load_axes('data/%s/%s'%(expr['path'],fname3))
        
        nlm_1[ff,:] = axes_to_nlm(axes1)
        nlm_2[ff,:] = axes_to_nlm(axes2)
        nlm_3[ff,:] = axes_to_nlm(axes3)
        
        # Deformation gradient
        fname = '%s-F.csv'%(fname0[:-4])
        df = pd.read_csv('data/%s/%s'%(expr['path'],fname), header=None, sep=',') 
        F[:,:,ff] = df.to_numpy()
                
    ### Save all correlations from the files of this experiment
    for ii, nlm in enumerate((nlm_1,nlm_2,nlm_3),1):
        fname = os.path.basename(expr['path'])[0:]
        fcorr = "drex-state-trajectories/%s-%i.p"%(fname, ii) 
        print('=== Saving states: %s ==='%(fcorr))
        pickle.dump([nlm,F], open(fcorr, "wb"))
               
