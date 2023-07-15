#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2023

"""
Discrete Directors Method (DDM) for initially isotropic vector bundle of crystallographic axes
"""
 
import numpy as np
import copy, sys, time, code # code.interact(local=locals())

class DDM():

    """
    Discrete Directors Method (DDM) solver
    """

    def __init__(self, iota=+1, N=1000, v0=None):
    
        self.iota = iota
        self.N = N
        self.v = np.zeros((N,3))        

        if v0 is None:
            # Initialize with randomly distributed vector bundle
            for ii in np.arange(N):
                phi = np.random.uniform(0, 2*np.pi)
                theta = np.arccos(1-2*np.random.uniform(0, 1))
                x = np.sin(theta) * np.cos(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(theta)
                self.v[ii,:] = [x,y,z]
        elif len(v0.flatten()) > 3: # more than 3 vectors? set the passed list of init vectors
            self.N = v0.shape[0]
            self.v = v0 # (num, xyz)
        else:
            for ii in np.arange(self.N): self.v[ii,:] = v0


    def evolve(self, L, dt):
        
        D = (L + L.T)/2 # L is velocity gradient
        W = (L - L.T)/2

        v0 = self.v.copy()
        v0sq = np.einsum('ni,nj->nij', v0,v0) 
        
        dvdt = np.einsum('ij,nj->ni', W, v0) # W.v
        asymprod = np.einsum('nij,jk->nik', v0sq, D) - np.einsum('ij,njk->nik', D, v0sq)
        dvdt += self.iota*np.einsum('nij,nj->ni', asymprod, v0) # iota*(v^2.D - D.v^2).v
        
        self.v = v0 + dt*dvdt
        
        # re-normalize
        norms = np.linalg.norm(self.v, axis=1)
        self.v = np.array([self.v[ii,:]/norms[ii] for ii in np.arange(self.N)])


    def get_eigenframe(self):
    
        a2_n = np.einsum('ni,nj->nij', self.v, self.v)
        a2 = np.mean(a2_n, axis=0)
        print(a2)
        eigvals, eigvecs = np.linalg.eig(a2)
        I = eigvals.argsort() 
        return eigvecs[:,I], eigvals[I] # last entry is largest eigenvalue pair
        
    
    def get_principal_dir_2D(self, outofplane='y'):

        eigvecs, eigvals = self.get_eigenframe()
#        print('--------')
#        print(eigvals)
#        print(eigvecs)
        if outofplane == 'y': 
            e1 = eigvecs[[0,2],-1] # x,z components of large eigenvalue vector
            e1 /= np.linalg.norm(e1)
            return e1
        else: return None
