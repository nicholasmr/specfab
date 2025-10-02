# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-

"""
Vertical stack of horizontally homogeneous layers of anisotropic polycrystalline ice
"""

import numpy as np

from .layer import *
from .layer_PPJ import *

class LayerStack:
    
    def __init__(self, nlm, z, beta=[0,], f=100e6, \
                 epsa=[3.17-0.034], epsc=[3.17], sigma=[0e-5], mu=1, \
                 modeltype='GTM', VERBOSE=0): 

        self.modeltype = modeltype # GTM (General Transfer Matrix) or FP (Fujita--Paren)
        if self.modeltype not in ['GTM','FP']: raise ValueError('Argument "modeltype" must be either "GTM" or "FP"')

        self.beta = beta # horizontal rotations (radians)
        self.f    = f    # radar frequency
        
        # Dimensions
        self.N_layers,_ = nlm.shape # number of layers 
        self.N_frames   = len(self.beta) # number of rotated frames
        self.layerlist  = np.arange(self.N_layers)
        self.framelist  = np.arange(self.N_frames)
        
        self.z  = np.hstack(([np.abs(z[1]-z[0])],z)) # Add isotropic surface layer with thickness equal to first subsurface layer 
        self.dz = np.abs(np.diff(self.z))
     
        # Layer parameters
        self.mu    = mu  
        self.epsa  = np.full((self.N_layers), epsa[0])  if len(epsa)==1  else epsa 
        self.epsc  = np.full((self.N_layers), epsc[0])  if len(epsc)==1  else epsc 
        self.sigma = np.full((self.N_layers), sigma[0]) if len(sigma)==1 else sigma 
        
        # Save spectral fabric profile
        self.nlm = nlm[:, :6] # Save only L=2 truncation (model does not depend on higher-order modes)
        self.lm  = np.array([(0,0), (2,-2),(2,-1),(2,0),(2,1),(2,2)]).T
        
        ### Errors, warnings and verbosity
        
        if len(self.z) != 1+len(self.nlm): raise ValueError('len(z) (=%i) != 1+len(nlm) (=%i)'%(len(self.z), 1+len(self.nlm)))
        
        if np.any(nlm[0,1:]): raise ValueError('n_2^m of top layer must vanish (surface layer must be isotropic)')
        
        self.VERBOSE = VERBOSE
        if self.VERBOSE: print('Initialized "%s" layer stack: frames x layers = %i x %i'%(self.modeltype, self.N_frames,self.N_layers))

        zeta = 0 # Dimensionless x-component of incident wave vector
        self.frm_layers = [self.get_layerstack(b, self.f, zeta) for b in self.beta]


    def get_layerstack(self, beta, f, zeta):
        return [Layer(self.nlm[nn], beta, self.dz[nn], f, zeta, self.epsa[nn], self.epsc[nn], self.sigma[nn], self.mu, modeltype=self.modeltype) for nn in self.layerlist] 


    def get_propconst(self):
        # Dimensionless propagation constants: k_s = omega/c q_s 
        return np.array([self.frm_layers[0][nn].ks[:] for nn in self.layerlist[:]], dtype=np.complex128) 
   
 
    ### RADAR RETURNS
    
    def get_returns(self, E0=1, Tx_pol=([1,0],), nn=None):
    
        ### Construct URD matrix for total two-way propagation
        
        URD = self.get_URD(self.frm_layers)
        
        ### Calculate wave returns

        if self.modeltype=='FP': E0 /= 13 # rescale FP transmission strength to match GTM (needs to be understood why, but does not affect polarimitric results)
    
        # Tx_pol:  polarization of transmitted wave H (x'' dir), V (y'' dir)
#        Tx_pol = ([1,0],)
        E_Tx = np.array([np.array([E0*p[0], E0*p[1]], dtype=np.complex128) for p in Tx_pol]) # Init downard propergating wave in layer 0
        E_Rx = self.get_Rx_all(E_Tx, URD, nn=nn) # (frame, layer, Tx/Rx pair)
        
        E_HH = E_Rx[:,:,0] # Tx,Rx = H,H
        E_HV = E_Rx[:,:,1] # Tx,Rx = H,V
        
        E_HH_abs = np.abs(E_HH)
        E_HV_abs = np.abs(E_HV)
        
        P_HH  = self.dB(E_HH_abs)
        Pm_HH = self.dB(E_HH_abs.mean(axis=0, keepdims=True))
        dP_HH = P_HH - Pm_HH
        
        P_HV  = self.dB(E_HV_abs)
        Pm_HV = self.dB(E_HV_abs.mean(axis=0, keepdims=True))
        dP_HV = P_HV - Pm_HV

        h,v = 0,1 # p (H), s (V) components
        I = self.layerlist[:-1] if nn==None else np.array([nn]) # Consider only the requested layer (nn) ?
        numer = np.multiply(       URD[:,I,h,h],  np.conjugate(URD[:,I,v,v]))
        denom = np.multiply(np.abs(URD[:,I,h,h]),       np.abs(URD[:,I,v,v]))
        c_HHVV = np.divide(numer, denom)
        
        # phase offset: c_HHVV is defined only up to an arbitrary phase shift <=> exp[i*(phi1+const.)]*exp[-i*(phi2+const.)] = exp[i*(phi1-phi2)]
        if self.modeltype=='GTM': 
#            c_HHVV *=  np.exp(-1j*np.pi)
            c_HHVV = -np.abs(c_HHVV)*np.exp(-1j*np.angle(c_HHVV)) # @TODO go through GTM matrices to determine where the phase difference compared to FP occurs.
        
        returns = (np.squeeze(Pm_HH), np.squeeze(Pm_HV), \
                   np.squeeze(dP_HH.T), np.squeeze(dP_HV.T), \
                   np.squeeze(c_HHVV.T), )

        return returns
        
    def dB(self, amp): return 20 * np.log10(amp) 

    def get_Rx_all(self, E_Tx_list, URD, nn=None):

        Nl = self.N_layers-1 if nn is None else 1
        E_Rx = np.zeros((self.N_frames, Nl, 2*len(E_Tx_list)), np.complex128) 

        for jj, E_Tx in enumerate(E_Tx_list): # Transmitted polarization
            I = [0+2*jj, 1+2*jj]        
            if nn is None: E_Rx[:,:,I] = np.einsum('fnij,j->fni', URD[:,:,:,:],  E_Tx) 
            else:          E_Rx[:,0,I] = np.einsum('fij,j->fi',   URD[:,nn,:,:], E_Tx) 
        
        return E_Rx
 
    """
    UNDER THE HOOD
    """
    
    ### SYSTEM MATRICES

    def get_T_R(self, Aprev,Aprev_inv, Anext,Anext_inv):
        
        # Calculates the 2x2 matrices for reflection (R), downward transmission (T), and upward transmission (Trev)
        
        ### Downward propergating wave (modified Jeannin (2019) GTM code)
       
        Delta1234 = np.array([[1,0,0,0],
                              [0,0,1,0],
                              [0,1,0,0],
                              [0,0,0,1]])
    
        Delta1234_inv = exact_inv(Delta1234)  
        
        Gamma     = np.einsum('fnij,fnjk->fnik', Aprev_inv, Anext) 
        GammaStar = np.einsum('ij,fnjk,kl->fnil', Delta1234_inv, Gamma, Delta1234) 
        
        # Common denominator for all coefficients
        Denom = np.multiply(GammaStar[:,:,0,0],GammaStar[:,:,2,2]) - np.multiply(GammaStar[:,:,0,2],GammaStar[:,:,2,0])
        
        # Reflection coefficients
        rpp = np.multiply(GammaStar[:,:,1,0],GammaStar[:,:,2,2]) - np.multiply(GammaStar[:,:,1,2],GammaStar[:,:,2,0])
        rss = np.multiply(GammaStar[:,:,0,0],GammaStar[:,:,3,2]) - np.multiply(GammaStar[:,:,3,0],GammaStar[:,:,0,2])
        rps = np.multiply(GammaStar[:,:,3,0],GammaStar[:,:,2,2]) - np.multiply(GammaStar[:,:,3,2],GammaStar[:,:,2,0])
        rsp = np.multiply(GammaStar[:,:,0,0],GammaStar[:,:,1,2]) - np.multiply(GammaStar[:,:,1,0],GammaStar[:,:,0,2])

        rpp = np.nan_to_num(np.divide(rpp,Denom))
        rss = np.nan_to_num(np.divide(rss,Denom))
        rps = np.nan_to_num(np.divide(rps,Denom))        
        rsp = np.nan_to_num(np.divide(rsp,Denom))
        
        R = np.array([[rpp,rsp], [rps,rss]], dtype=np.complex128) 

        # Transmission coefficients
        tpp = np.nan_to_num(+np.divide(GammaStar[:,:,2,2],Denom))
        tss = np.nan_to_num(+np.divide(GammaStar[:,:,0,0],Denom))
        tps = np.nan_to_num(-np.divide(GammaStar[:,:,2,0],Denom))
        tsp = np.nan_to_num(-np.divide(GammaStar[:,:,0,2],Denom))

        T = np.array([[tpp, tsp], [tps, tss]], dtype=np.complex128)
        
        ### Upward propergating wave
        
        Gamma_rev = np.einsum('fnij,fnjk->fnik', Anext_inv, Aprev)

        # Common denominator for all coefficients
        Denom_rev = np.multiply(Gamma_rev[:,:,2,2],Gamma_rev[:,:,3,3])  - np.multiply(Gamma_rev[:,:,2,3],Gamma_rev[:,:,3,2]) 

        # Transmission coefficients
        tpp_rev = np.nan_to_num(+np.divide(Gamma_rev[:,:,3,3],Denom_rev))
        tss_rev = np.nan_to_num(+np.divide(Gamma_rev[:,:,2,2],Denom_rev))
        tps_rev = np.nan_to_num(-np.divide(Gamma_rev[:,:,3,2],Denom_rev))
        tsp_rev = np.nan_to_num(-np.divide(Gamma_rev[:,:,2,3],Denom_rev))

        Trev = np.array([[tpp_rev, tsp_rev], [tps_rev, tss_rev]], dtype=np.complex128)
        
        # Return re-ordered tensors for convention used in code (frame,layer, ...)
        return np.einsum('ijfn->fnij',R), np.einsum('ijfn->fnij',T), np.einsum('ijfn->fnij',Trev)

        
    def get_URD(self, frm_layers):

        D = np.zeros((self.N_frames,self.N_layers, 2,2), dtype=np.complex128) # Downward (fwd) propagation
        U = np.zeros((self.N_frames,self.N_layers, 2,2), dtype=np.complex128) # Upward (rev) propagation
        D[:,0,:,:] = np.eye(2) # identity for surface layer (zero propagation path)
        U[:,0,:,:] = np.eye(2) # identity for surface layer (zero propagation path)
        
        laylist = self.layerlist[:-1] # Reduced layer list: can not calculate reflection from last layer since the interface matrix depends on the perm/fabric of the next (unspecified) layer

        ### Construct layer/interface-wise transfer matrices

        Mfwd_next = np.array([ [frm_layers[ff][nn+1].Mfwd for nn in laylist] for ff in self.framelist], dtype=np.complex128)        
        Mrev_next = np.array([ [frm_layers[ff][nn+1].Mrev for nn in laylist] for ff in self.framelist], dtype=np.complex128)

        if self.modeltype == 'GTM':
            
            A         = np.array([ [frm_layers[ff][nn  ].Ai     for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            A_inv     = np.array([ [frm_layers[ff][nn  ].Ai_inv for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            Anext     = np.array([ [frm_layers[ff][nn+1].Ai     for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            Anext_inv = np.array([ [frm_layers[ff][nn+1].Ai_inv for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            
            R, T, Trev = self.get_T_R(A,A_inv, Anext,Anext_inv)
            
        if self.modeltype == 'FP':
    
            R    = np.zeros((self.N_frames,len(laylist), 2,2), dtype=np.complex128)
            T    = np.zeros((self.N_frames,len(laylist), 2,2), dtype=np.complex128)
            Trev = np.zeros((self.N_frames,len(laylist), 2,2), dtype=np.complex128)
            
            # Transmission matrices = identity matrices for all layers and frames
            T[:,:,0,0] = 1 # x,x entry
            T[:,:,1,1] = 1 # y,y entry
            Trev = T.copy()
            
            # Reflection matrices are based on Paren (1981) scattering
            ff = 0 # any beta frame will do (eigenvalues are not changed by rotation the frame of reference)
            R0flat = np.array([[ \
                            (frm_layers[ff][nn].ai_horiz[ii] - frm_layers[ff][nn+1].ai_horiz[ii])*(self.epsc[nn]-self.epsa[nn]) \
                      for ii in (0,1)] for nn in laylist], dtype=np.complex128)
            #R0flat *= 1/0.03404121647756628 # Cal. against GTM ?
            R0 = np.zeros((len(laylist), 2,2), dtype=np.complex128) # proper shape, not flattened
            R0[:,0,0] = +R0flat[:,0]
            R0[:,1,1] = +R0flat[:,1]
            # ...construct reflection matrices for each frame by a simple rotation around the z-axis or the reference frame
            for ff,b in enumerate(self.beta):
                for nn in self.layerlist[:-1]:
                    R[ff,nn,:,:] = frm_layers[ff][nn+1].eig_to_meas_frame(R0[nn,:,:])
        
        ### Construct total upward and downward transfer matrices for scattering of each internal interface.
        
        D_mul = np.einsum('fnij,fnjk->fnik', Mfwd_next, T)
        U_mul = np.einsum('fnij,fnjk->fnik', Trev, Mrev_next)
        
        for nn in self.layerlist[:-1]:
            D[:,nn+1,:,:] = np.einsum('fij,fjk->fik', D_mul[:,nn,:,:],     D[:,nn,:,:]) # D_{n}*E_{0}^{-,Tx} = E_{n}^{-,downward}
            U[:,nn+1,:,:] = np.einsum('fij,fjk->fik',     U[:,nn,:,:], U_mul[:,nn,:,:]) # U_{n}*E_{n}^{-,upward} = E_{0,Rx}^{-}

        ### Total downwards-reflected-upwards transformation (self.N_frames,self.N_layers-1, 2,2)
        
        URD = np.einsum('fnij,fnjk->fnik',U[:,:-1,:,:],np.einsum('fnij,fnjk->fnik',R,D[:,:-1,:,:])) 
        
        return URD
    
    def get_U(self, frm_layers):
    
        U = np.zeros((self.N_frames,self.N_layers, 2,2), dtype=np.complex128) # Upward (rev) propagation
        U[:,0,:,:] = np.eye(2) # identity for surface layer (zero propagation path)
        U[:,-1,:,:] = np.eye(2) # identity for bottom (tranmission) layer (zero propagation path)
        
        laylist = self.layerlist[:-1] # Reduced layer list: can not calculate reflection from last layer since the interface matrix depends on the perm/fabric of the next (unspecified) layer

        ### Construct layer/interface-wise transfer matrices

        Mrev_next = np.array([ [frm_layers[ff][nn+1].Mrev for nn in laylist] for ff in self.framelist], dtype=np.complex128)

        if self.modeltype == 'GTM':
            A         = np.array([ [frm_layers[ff][nn  ].Ai     for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            A_inv     = np.array([ [frm_layers[ff][nn  ].Ai_inv for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            Anext     = np.array([ [frm_layers[ff][nn+1].Ai     for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            Anext_inv = np.array([ [frm_layers[ff][nn+1].Ai_inv for nn in laylist] for ff in self.framelist], dtype=np.complex128)
            R, T, Trev = self.get_T_R(A,A_inv, Anext,Anext_inv)
        else:
            raise ValueError('Only GTM model type supported')
            
        ### Construct total upward and downward transfer matrices for scattering of each internal interface.
        
        U_mul = np.einsum('fnij,fnjk->fnik', Trev, Mrev_next)
        for nn in self.layerlist[:-1]:
            U[:,nn+1,:,:] = np.einsum('fij,fjk->fik', U[:,nn,:,:], U_mul[:,nn,:,:]) # U_{n}*E_{n}^{-,upward} = E_{0,Rx}^{-}

        return U[:,:-1,:,:]
    
