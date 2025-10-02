# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-

"""
Wraps a spectral fabric formulation around the Passler and Paarmann (2017,2019) layer class implemented by M. Jeannin (2019)
"""

from ..specfabpy import specfabpy as sf__

import numpy as np
import numpy.linalg as lag

from .layer_PPJ import *

class Layer(Layer_PPJ):
    
    def __init__(self, nlm, beta, thickness, f, zeta, epsa, epsc, sigma, mu, modeltype='GTM'): 
        
        super().__init__()
        
        self.modeltype = modeltype # GTM (General Transfer Matrix) or FP (Fujita--Paren)
        if self.modeltype not in ['GTM','FP']: raise ValueError('Argument "modeltype" must be either "GTM" or "FP"')
        
        self.thick  = thickness # Layer thickness (d)
        self.f      = f         # Plane wave frequency 
        self.zeta   = zeta      # x-component of wave vector (dimensionless)
        
        self.epsa  = epsa  # Complex permitivity in plane perdendicular to c-axis
        self.epsc  = epsc  # Complex permitivity parallel to c-axis
        self.sigma = sigma # Single crystal isotropic conductivity
        self.mu    = mu    # Relative permeability
        
        # Short hands
        self.omega = 2*np.pi*self.f     # Angular frequency of wave
        self.q2k   = self.omega/c_const # Convert from q to k 
                
        self.setup(nlm, beta) 

        
    def setup(self, nlm, beta):
        
        # Given the orientation fabric (nlm), all layer matrices are calculated by this method.
        # Call this method if you wish to update the layer with a new orientation fabric and re-calculate all layer properties (matrices etc.).
        
        # Fabric frame 
        self.nlm_0   = nlm # Save copy
        self.ccavg_0 = sf__.a2(self.nlm_0) # <c^2> in fabric frame
        
        # Measurement frame
        self.beta = beta
        Qz = np.diag([np.exp(+1j*m*self.beta) for m in [0, -2,-1,0,+1,+2]]) # Rotation matrix for rotations about the vertical z-axis.
        self.nlm = np.matmul(Qz, self.nlm_0) 
        self.ccavg = sf__.a2(self.nlm) # <c^2> in measurement frame
        
        # Calculate permitivity tensor and (transfer) matrices
        self.epsilon = self.get_epsavg(self.ccavg)  
        self.set_matrices() # Requires self.epsilon be set
        

    def get_epsavg(self, ccavg):
    
        # Bulk (grain averaged) dielectric permittivity tensor from second-order structure tensor, <c^2> 
        
        monopole   = (2*self.epsa+self.epsc)/3 * np.eye(3)
        quadrupole = (self.epsc-self.epsa) * (ccavg - np.eye(3)/3)
        epsavg     = np.real(monopole+quadrupole) # Should be real; numerical errors might cause small non-zero imaginary parts.
        
        epsdprime = self.sigma/(2*np.pi*self.f*eps0) 
        epsavg = epsavg - 1j*epsdprime * np.eye(3) # Add isotropic imaginary part (due to acidity)
        
        return epsavg
    
    
    def set_matrices(self):
        
        # Sets matrices etc. needed to construct layer stack (LayerStack class)
        
        ### Horizontal eigen frame 
        
        self.ai_horiz, self.Qeig = lag.eig(self.ccavg_0[:2,:2])
        self.Qeiginv = self.Qeig.T # or lag.inv(self.Qeig)
        self.eigang = np.arctan2(self.Qeig[1,0],self.Qeig[0,0]) # Horizontal rotation of fabric to get from measurement frame to eigen frame
        
        ### Setup layer matrices
        
        # General Transfer Matrix (Rathmann et al., 2022)
        if self.modeltype == 'GTM': 
        
            # Call parent GTM routines
            self.calculate_matrices(self.zeta)
            self.calculate_q()
            self.ks = self.q2k*self.qs # Propagation constants (not used in GTM, used for comparing with FP model)
            self.calculate_gamma(self.zeta)
            self.calculate_transfer_matrix(self.f, self.zeta)  
            
            # Calculate propagation and auxiliary matrices
            self.Ai_inv = exact_inv(self.Ai.copy())
            self.Ki_inv = exact_inv(self.Ki.copy()) # See Layer_PPJ: diag(exp(-1j*(omega/c)*qs[:]*d)) = diag(exp(-1j*ks[:]*d)
            self.Mfwd = self.Ki_inv[0:2,0:2] # downwards: +1j*k*d exponent
            self.Mrev = self.Ki[2:,2:]       # upwards:   -1j*k*d exponent (Yeh/Passler sign convention)
        
        # Fujita--Paren
        if self.modeltype == 'FP': 
            
            # Bulk dielectric permittivity in x,y directions
            # Diagonal components of eqn. (4) in Fujita et al. (2006)
            epsxy = self.epsa + self.ai_horiz*(self.epsc-self.epsa) # Horizontal components when <cc> is in eigen frame (see eqn. (2) in Rathmann et al., 2021)
            
            # Propagation constants kx and ky in Fujita et al. (2006), eqn. (7a,7b)
            kxy = np.sqrt(epsxy*(self.q2k)**2 + 1j*mu0*self.sigma*self.omega) 
#            kxy = -np.conj(kxy) # debug: for changing polarity of coherence phase c_HHVV.
            (kx,ky) = kxy # Decompose for clarity
            self.ks = +np.array([kx,ky, kx,ky], dtype=np.complex128)
            self.qs = 1/self.q2k*self.ks
            
            # Calculate propagation matrices for layer in eigen frame
            # "Transmission" matrix, eqn. (5) in Fujita et al. (2006)
            self.Mfwd = np.diag( np.exp(+1j*self.ks[[0,1]]*self.thick) )
            self.Mrev = np.diag( np.exp(+1j*self.ks[[2,3]]*self.thick) )
            
            # Rotate eigen frame back to the measurement frame
            self.Mfwd = self.eig_to_meas_frame(self.Mfwd)
            self.Mrev = self.eig_to_meas_frame(self.Mrev)


    def eig_to_meas_frame(self, mat_eig_frame):
        
        ang = -self.beta + self.eigang # Fabric eigen frame to measurement frame
        c, s = np.cos(ang), np.sin(ang)
        Qz = np.array(((c, -s), (s, c))) # rotation matrix
        return np.matmul(Qz, np.matmul(mat_eig_frame, Qz.T)) # = mat_meas_frame


