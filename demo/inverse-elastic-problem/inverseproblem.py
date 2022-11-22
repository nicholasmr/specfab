# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Inversion class that can
    1) Infer monocrystal elastic parameters ("g" vector) given an ODF (and phase-velocity measurements)
    2) Infer harmonic coefficients of an ODF ("nlm" vector) given monocrystal elastic parameters (and phase-velocity measurements)
"""

import copy, sys, code # code.interact(local=locals())
import math
import numpy as np 
import numpy.linalg as linalg
from scipy.optimize import minimize
sys.path.insert(0, '..')
from specfabpy import specfabpy as sf 

"""
Inverse problem class
"""

def f_J(misfitvec, beta, N): 
    # Cost-function J (data--model velocity misfit measure)
    return beta[0]/N*np.sum(misfitvec[0]**2) + beta[1]/N*np.sum(misfitvec[1]**2) + beta[2]/N*np.sum(misfitvec[2]**2) # Terms are the P, S1, and S2 misfits, respectively

def x_to_nlm(x):
    # Vector of independent nlm coefficients to infer, collapsed to 
    #   an array of real numbers for the inversion routine
    n20  = x[0]
    n21p = x[1] + 1j*x[2]
    n22p = x[3] + 1j*x[4]
    n40  = x[5] 
    n41p = x[6] + 1j*x[7]
    n42p = x[8] + 1j*x[9]
    n43p = x[10] + 1j*x[11]
    n44p = x[12] + 1j*x[13]
    return np.array([1/np.sqrt(4*np.pi), \
                                            +np.conj(n22p), -np.conj(n21p), n20, n21p, n22p, \
            +np.conj(n44p), -np.conj(n43p), +np.conj(n42p), -np.conj(n41p), n40, n41p, n42p, n43p, n44p], dtype=np.complex128)

class InverseProblem:
   
    def __init__(self, rho=917, gtol0=1e-3, verbose=False):
        
        self.lm, self.nlm_len = sf.init(4) # L=4 suffices here
        self.rho = rho # Mass density

        ### Inversion numerics
        self.method = 'CG'
        self.gtol0  = gtol0  # tolerence for p
        self.gtol1  = self.gtol0 # tolerence for nlm
        self.verbose = verbose

    def infer_g(self, observations, alpha, beta, g0):
    
        """
        Infer the grain elastic parameters (g vector)
        """

        # Cost function        
        def J_g(g, *args):
            vi_obs,nlm_obs,theta,phi, alpha,rho, beta,scalevec = args
            (vP, vS1, vS2) = get_vjmap(nlm_obs, alpha, scalevec*g, rho, theta, phi)
            (vP_obs,vS1_obs,vS2_obs) = vi_obs # unpack
            return f_J((vP-vP_obs, vS1-vS1_obs, vS2-vS2_obs), beta, len(phi)) 

        (vi_obs,nlm_obs,theta,phi) = observations
        scale0 = 1.0e10
        scalevec = np.array([scale0,scale0, 1,1,1]) # scale "p" during cost-function minimization to not encounter numerical accuracy errors (lam and mu are ~1e9, whereas Ei are ~ 1e0)        
        res = minimize(J_g, g0/scalevec, (vi_obs,nlm_obs,theta,phi, alpha, self.rho, beta, scalevec), method=self.method, options={'disp':self.verbose, 'gtol': self.gtol0})
        g = scalevec*res.x
        vi_g = get_vjmap(nlm_obs, alpha, g, self.rho, theta, phi) # predicted velocities using best fit params
        dvi_g = get_dvi(vi_g)
        return (g, vi_g, dvi_g)


    def infer_nlm(self, observations, alpha, g, beta, eta, use_angular_anomalies=False):
    
        """ 
        Infer the ODF (\hat{psi})
        """

        # Cost function
        def J_nlm(x, *args):
            vi_obs,theta,phi, alpha,g,rho, beta,eta, use_angular_anomalies = args    
            nlm = x_to_nlm(x[:-2])
            s1,s2 = x[-2],x[-1] # slack variables
            vi = get_vjmap(nlm, alpha,g,rho, theta,phi)
            misfit_vector = get_dvi(vi) - get_dvi(vi_obs) if use_angular_anomalies else vi - vi_obs # dvi = delta vi = vi - <vi>
            J_V = f_J(misfit_vector, beta, len(phi))
                        
            ## Inequality constraints added with slack variables
            lam1 = np.amin(linalg.eigvals(sf.a2(nlm))) # smallest <c^2> eigenvalue
            _,_,_,_,_,_, Lami = sf.a4_eigentensors(nlm)
            Lam1 = np.amin(Lami) # smallest <c^4> eigenvalue 
            J_s = -eta[0]*(lam1-s1**2) - eta[1]*(Lam1-s2**2)
            
            ## Additional regularization
            J_r = 0 # No additional regularization
            #J_r = 0.5e2*np.sum(np.abs(nlm)) # Tikhonov regularization (solutions with smaller vector norms are preferred). 

            ## Construct total cost
            J = J_V + J_s + J_r
            return np.log(J) # seems to work slightly better taking log(J) compared to just J

        (vi_obs, theta,phi) = observations
        x0 = np.zeros((3+2 + 5+4 + 2)) # init guess: nlm=0 => isotropy
        res = minimize(J_nlm, x0, (vi_obs,theta,phi, alpha,g,self.rho, beta,eta, use_angular_anomalies), \
                            method=self.method, options={'disp':self.verbose, 'gtol': self.gtol1,})
        nlm_infr = x_to_nlm(res.x[:-2]) # two last values are the slack variables
        vi_infr = get_vjmap(nlm_infr, alpha, g, self.rho, theta, phi) # predicted velocities using best fit ODF
        dvi_infr = get_dvi(vi_infr)
        return (nlm_infr, vi_infr, dvi_infr)

    
"""
Shared structure and routines
"""

lm_L4 = np.array([(0,0), (2,-2),(2,-1),(2,0),(2,1),(2,2), (4,-4),(4,-3),(4,-2),(4,-1),(4,0),(4,+1),(4,+2),(4,+3),(4,+4)]).T 

def get_vjmap(nlm, alpha, g, rho, theta, phi, freq=10):
    (lam,mu,Elam,Emu,Egam) = g
    omega = 2*np.pi * freq # 2*pi * freq, but freq does not matter for results!
    vj_flat = sf.Vi_elastic(nlm, alpha, lam,mu,Elam,Emu,Egam, omega, rho, theta,phi)
    vS1_true = vj_flat[0,:]
    vS2_true = vj_flat[1,:]
    vP_true  = vj_flat[2,:]
    return np.array([vP_true, vS1_true, vS2_true])

def get_dvi(vi):
    # dvi = delta vi = vi - <vi>
    (vP, vS1, vS2) = vi
    dvP  = vP  - np.mean(vP)
    dvS1 = vS1 - np.mean(vS1)
    dvS2 = vS2 - np.mean(vS2)
    return np.array([dvP, dvS1, dvS2])
    
def cart2sph(v, deg=False, colat=False):
    x,y,z = v[0],v[1],v[2]
    theta = math.atan2(math.sqrt(x**2 + y**2), z)
    phi = math.atan2(y, x) #if x >= 0 else math.atan2(y, x) + math.pi
    if phi < 0: phi = 2*np.pi + phi
    if colat: theta = np.pi/2 - theta # 
    if deg: theta, phi = np.rad2deg(theta), np.rad2deg(phi)
    return theta, phi
    
