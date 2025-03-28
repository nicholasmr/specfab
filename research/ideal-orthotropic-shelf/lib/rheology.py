#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024

import numpy as np 
from .constants import *
from specfabpy.fenics.rheology import *

class Rheology(ASSA):

    """
    Rheology
    """
    
    def __init__(self, **kwargs):
    
        ### Save arguments
        self.__dict__.update(kwargs)
        
        ### Init
        super().__init__(**kwargs)
