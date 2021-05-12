#!/usr/bin/python
# N. M. Rathmann, 2019

import numpy as np

#-------------------
# PURE SHEAR (PS)
#-------------------

class PureShear():

    def __init__(self, t_e, r, vertdir='z'): 

        self.t_e = float(t_e) # e-folding time scale

        if vertdir=='z': self.Fpow = [(1+r)/2., (1-r)/2., -1]
        if vertdir=='x': self.Fpow = [-1, (1+r)/2., (1-r)/2.]

    def lam(self, t):    return np.exp(t/self.t_e) # lambda(t)
    def F(self, t):      return np.diag(np.power(self.lam(t),self.Fpow))
    def strain(self, t): return 0.5*( self.F(t) + np.transpose(self.F(t)) ) - np.diag([1,1,1])

    def epszz2time(self,epszz): return -self.t_e*np.log(epszz+1) # time it takes to reach "eps_zz" strain with t_e char. timescale

    # Note that F is constructed such that W and eps are time-independant.
    def W(self):   return np.diag([0,0,0])
    def eps(self): return 1./self.t_e * np.diag(self.Fpow)


#-------------------
# SIMPLE SHEAR (SS) 
#-------------------

class SimpleShear():

    def __init__(self, kappa, shearplane='xz'): 

        self.kappa0 = kappa    
        self.kappa = [0,0]

        if shearplane=='xz': self.kappa[0] = self.kappa0
        if shearplane=='xy': self.kappa[1] = self.kappa0

    def time2shearangle(self, t): return np.arctan(self.kappa0*t)
    
    def shearangle2time(self, ang): return np.tan(ang)/self.kappa0
    
    # Note that F is constructed such that W and eps are time-independant.
    def W(self):   return 0.5*np.matrix([[0,+self.kappa[1],+self.kappa[0]], [-self.kappa[1],0,0], [-self.kappa[0],0,0]]);
    def eps(self): return 0.5*np.matrix([[0,+self.kappa[1],+self.kappa[0]], [+self.kappa[1],0,0], [+self.kappa[0],0,0]]);

#-------------------
# RIGID ROTATION (RR) 
#-------------------

class RigidRotation():

    def __init__(self, omega, rotax='z'): 

        self.omega0 = omega
        self.omega = [0,0,0]

        if rotax=='x': self.omega[0] = self.omega0
        if rotax=='y': self.omega[1] = self.omega0
        if rotax=='z': self.omega[2] = self.omega0

    def beta(self, t): return self.omega0*t

    # Note that F is constructed such that W and eps are time-independant.
    def W(self):   return np.matrix([[0,-self.omega[2],self.omega[1]], [self.omega[2],0,-self.omega[0]], [-self.omega[1],self.omega[0],0]]);
    def eps(self): return np.diag([0,0,0])


