#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Default experiment configuration
"""

class Config():


    ### Mesh

    modelplane = 'xz'
    geom       = 'rectangular' # supported: "rectangular"
    h0profile  = 'triangular'
    H0         = 500  # shelf height at inflow
    evolve     = True # for debugging

    def kwargs_mesh(self): 
        return dict(modelplane=self.modelplane, geom=self.geom, h0profile=self.h0profile, H0=self.H0, evolve=self.evolve,)


    ### Bulk flow law

    rheology = 'Orthotropic' # supported: "Isotropic", "Orthotropic"
    n        = 3             # flow exponent
    T        = -15 + 273.15  # isothermal temperature (deg. K)

    def kwargs_rheo(self): 
        return dict(rheology=self.rheology, n=self.n, T=self.T, )


    ### Fabric evolution and homogenization

    fabdyn        = 'LROT'   # Crystal process modelled (LROT or DDRX)    
    L             = 10       # spectral truncation
    nu_realspace  = 1e-5     # real-space regularization 
    nu_multiplier = 0.9      # S^2 regularization multiplier
    
    alpha         = 0.455    # Taylor--Sachs homogenization weight
    Eij_grain     = (1, 1e3) # (Ecc, Eca) grain enhancements
    n_grain       = 1        # grain power-law exponent (only n_grain=1 supported)
    DDRX_symframe = 4        # a4 eigentensor number (0 to 5) to use for estimating mi for DDRX fabrics
    
    def kwargs_fabric(self): 
        return dict(fabdyn=self.fabdyn, L=self.L, nu_realspace=self.nu_realspace, nu_multiplier=self.nu_multiplier, \
                    alpha=self.alpha, Eij_grain=self.Eij_grain, n_grain=self.n_grain, \
                    modelplane=self.modelplane, DDRX_symframe=self.DDRX_symframe, ) 

    ### All flowmodel kwargs 

    def kwargs_flowmodel(self): 
        return dict( kwargs_rheo=self.kwargs_rheo(), kwargs_fabric=self.kwargs_fabric(), kwargs_mesh=self.kwargs_mesh(), )

    ### Time stepping
    
    nt = 0 # Starting time step (used to resume integration if non-zero)
    Nt = 500 +1
    tt_savemod = 5 # number of steps between plotting diagnostics

    ### File paths 
    
    def path_output(self):      return 'experiments/%s-dyn=%s-L=%i'%(self.rheology, self.fabdyn, self.L)
    def path_diagnostics(self): return './%s/diagnostics'%(self.path_output())
    def path_statedump(self):   return './%s/statedump'%(self.path_output())

    def fname_state(self, tt):  return '%s/%04d.h5'%(self.path_statedump(), tt)
    def fname_diagn(self, tt):  return '%s/diagnostic%04i.png'%(self.path_diagnostics(), tt)

    ### Plot ODF at these points in diagnostic plots
    
    odf_x = [] 
    odf_y = []

