#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Experiment configuration file for Rathmann et al. (2025)
"""

from lib.defaultconfig import *

class R25config(Config):

    odf_x = [30e3, 40e3, 50e3, 60e3] 
    odf_y = [1e3, 5e3, 10e3, 14e3]

    def path_output(self): return 'experiments/%s-dyn=%s-L=%i'%(self.rheology, self.fabdyn, self.L)
