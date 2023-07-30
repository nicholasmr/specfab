# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Wrapper for experimental data sets. 
Provides angular phase-velocity measurements and harmonic expansion coefficients of ODF.
"""

import numpy as np 
import copy, sys, code # code.interact(local=locals())
import pandas as pd
import quaternion as qt
from scipy.optimize import minimize

from Cij import * # elastic constants
from inverseproblem import get_vi_map, f_J
sys.path.insert(0, '..')
from specfabpy import specfabpy as sf 

from plottools import cart2sph

class Lutz_etal_2022:

    PATH_VELOCITIES = 'data/Lutz_etal_2022' 
    PATH_CAXES      = 'data/Thomas_etal_2021' 
   
    def __init__(self, rho=917, verbose=False):
        
        self.lm, self.nlm_len = sf.init(4) # L=4 suffices here
       
        # Used for optimization problem in order to find ODF orientation relative to phase-velocity measurement frame
        self.rho   = rho # Mass density
        self.alpha = 0.5 # Hill average
        self.beta  = [1,1,1] # Misfit weights
        self.g = Cij_to_g(Cij['Bennett1968']) # = (lam,mu,gam, Elam,Emu,Egam)
        self.g = self.g[[0,1,3,4,5]] # Routines below depend only on independent parameters (recall that gamma=lambda+2*mu)

        self.method = 'CG'
        self.gtol0  = 1e-3  # tolerance for g
        self.verbose = verbose
       
    def get_velocities(self, exprnum):

        kwargs = {'delimiter':" ", 'comment':'%', 'names':["phi", "vi", "sig_vi"]}

        data = pd.read_csv('%s/%03i_vp_kinematic.txt'%(self.PATH_VELOCITIES,exprnum), **kwargs)
        vP_obs, vP_obs_sig = data["vi"].to_numpy(), data["sig_vi"].to_numpy()

        data = pd.read_csv('%s/%03i_vs1_kinematic.txt'%(self.PATH_VELOCITIES,exprnum), **kwargs)
        vS1_obs, vS1_obs_sig = data["vi"].to_numpy(), data["sig_vi"].to_numpy()
        
        data = pd.read_csv('%s/%03i_vs2_kinematic.txt'%(self.PATH_VELOCITIES,exprnum), **kwargs)
        vS2_obs, vS2_obs_sig = data["vi"].to_numpy(), data["sig_vi"].to_numpy()
        
        # Lutz et al.: 
        #  - phi is the sample orientation relative to macroscopic flow direction, i.e. phi=0 => parallel to shear margin
        #  - phi is defined as increasing in the clockwise direction (from +y-axis towards +x-axis), 
        #       whereas we adopt the usual polar coordinate definition (phi positive anti-clockwise)
#        phi = np.deg2rad(data["phi"].to_numpy()) # forward phi direction
        phi = 360-data["phi"].to_numpy() + 90 # reverse phi direction, zero along +y axis, increases towards +x axis
        phi[phi>360] -= 360 # fix overflow
        phi = np.deg2rad(phi) 
        I = np.argsort(phi) # re-sort data points for monotonically increasing phi 
        phi = phi[I]
        theta = np.pi/2 + 0*phi # measurements were made in the horizontal x--y plane

        vi_obs    = (vP_obs[I], vS1_obs[I], vS2_obs[I])
        visig_obs = (vP_obs_sig[I], vS1_obs_sig[I], vS2_obs_sig[I])
        
        return (theta,phi, vi_obs, visig_obs)
        
       
    def get_nlm(self, exprnum, weighted=False):
    
        # In horizontal frame that is consistent with phase velocity measurements (kinematic frame)
        
        (theta,phi, vi_obs, _) = self.get_velocities(exprnum)
        (vP_obs,vS1_obs,vS2_obs) = vi_obs # unpack
        
        def J_phi(deltaphi, *args):
            theta,phi, nlm_obs, g,rho,alpha,beta = args
            (vP,vS1,vS2) = get_vi_map(sf.rotate_nlm(nlm_obs, 0,deltaphi), alpha, g,rho, theta, phi)
            return f_J((vP-vP_obs, vS1-vS1_obs, vS2-vS2_obs), beta, len(phi)) 
        
        (nlm_obs, qlat,qlon) = self.get_nlm_raw(exprnum, weighted=weighted, deltaphi=0) # unrotated ODF 
        deltaphi0 = np.deg2rad(-45*1) # initial guess 
        res = minimize(J_phi, deltaphi0, (theta,phi, nlm_obs, self.g,self.rho,self.alpha,self.beta), method=self.method, options={'disp':self.verbose, 'gtol': self.gtol0})
        deltaphi = res.x # best fit 
        if 0:
            print('d phi = %f'% np.rad2deg(deltaphi))
            deltaphi = 0 # for debugging, don't rotate ODFs but show in EBSD measurement frame
        returnvector = self.get_nlm_raw(exprnum, weighted=weighted, deltaphi=deltaphi)
        return returnvector # (nlm_obs, qlat,qlon)    


    def get_nlm_raw(self, exprnum, weighted=False, deltaphi=0):
    
        # In coordinate system of Thomas et al. (2021)
        
        fname = '%03d.ctf.csv'%(exprnum)
        df = pd.read_csv('%s/%s'%(self.PATH_CAXES, fname))
        qs_comps = df.to_numpy()[:,0:4]
        qs = qt.as_quat_array(qs_comps) # c-axes are the z-axis subject to these (quaternion) rotations
        
        R = qt.as_rotation_matrix(qs)
        v = np.array([np.matmul(R[ii,:,:], [0,0,1]) for ii in np.arange(len(R))])
        v = np.vstack((v, -v)) # add reflected vectors
        sphcoords = np.array([cart2sph(v[ii,:]) for ii in np.arange(len(v))])
        qcolat, qlon = sphcoords[:,0], sphcoords[:,1]
        qlat = np.pi/2 - qcolat     
        qlon += deltaphi        
        
        if weighted:
            area = df.to_numpy()[:,4] # grain area
            area = np.concatenate((area,area)) # add reflected
            area /= np.sum(area[:]) # normalized area (so that it sums to 1)
        else:
            area = df.to_numpy()[:,3] # just to get correct array size..
            area[:] = 1/len(qlon) # equal weights: 1/N

        # Construct a^4 and derive nlm using specfab        
        nlm = np.zeros((self.nlm_len), dtype=np.complex128) # psi expansion coefficients, n_l^m (psi_l^m in paper)
        caxes = np.array([ [np.cos(p)*np.sin(t), np.sin(p)*np.sin(t), np.cos(t)] for t, p in zip(qcolat,qlon) ])
        if weighted: a4 = np.array([ area[ii]*np.einsum('i,j,k,l',c,c,c,c) for ii,c in enumerate(caxes)]).sum(axis=0)
        else:        a4 = np.array([          np.einsum('i,j,k,l',c,c,c,c) for ii,c in enumerate(caxes)]).mean(axis=0)
        nlm[:sf.L4len] = sf.a4_to_nlm(a4)
        
        return (nlm, qlat,qlon)
    
