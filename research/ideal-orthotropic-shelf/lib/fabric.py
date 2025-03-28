#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024

import numpy as np 
from .constants import *
#from .defaultsettings import *

from specfabpy.fenics.ice import *

class Fabric(IceFabric):
    
    """
    Fabric evolution
    """

    def __init__(self, mesh, *args, **kwargs):

        ### Save arguments
        self.__dict__.update(kwargs)
    
        ### Set fields
        if self.fabdyn not in ['LROT','DDRX']: 
            raise ValueError('Fabric(): "fabdyn" must be either LROT or DDRX')
#        self.symframe = -1 if self.fabdyn=='LROT' else 4
        self.symframe = -1
               
        ### Init
        super().__init__(mesh.mesh, mesh.boundaries, symframe=self.symframe, **kwargs)
        self.set_isotropic_BCs([DOM_ID__IN,]) # isotropic CPO on inflow boundary

    def evolve(self, u, S, T, dt, *args, **kwargs):
        if self.fabdyn == 'LROT': kwargs_fabdyn = dict(iota=+1, Gamma0=None) 
        if self.fabdyn == 'DDRX': kwargs_fabdyn = dict(iota=+1, Gamma0=self.Gamma0_Lilien23_lab(u, T))
        super().evolve(u, S, dt, *args, **kwargs_fabdyn) # updates fabric state vector

