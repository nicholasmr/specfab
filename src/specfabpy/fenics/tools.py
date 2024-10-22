#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024

import copy, sys, time, code # code.interact(local=locals())
 
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from dolfin import *

from .. import constants as sfconst

def velgrid2velfuncspace(xv,yv, ux,uy, V=None, \
                          method='nearest', bounds_error=False, fill_value=None):
    
    """ 
    Contructs interpolants for populating FE velocity function space.
    NOTE: Assumes 2D problem.
    """

    # Interpolants
    points = (xv,yv)
    kwargs = dict(method=method, bounds_error=bounds_error, fill_value=fill_value)
    uxi = RegularGridInterpolator(points, ux, **kwargs)
    uyi = RegularGridInterpolator(points, uy, **kwargs)

    # Construct velocity function space
    class VelExpr(UserExpression):
        def __init__(self, uxi, uyi, **kwargs):
            super().__init__(**kwargs) 
            self.uxi = uxi
            self.uyi = uyi
        def eval(self, value, x):
            value[0] = self.uxi([x[0],x[1]])[0]
            value[1] = self.uyi([x[0],x[1]])[0]
        def value_shape(self):
            return (2,)
             
    uexpr = VelExpr(uxi, uyi)
    u = uexpr if V is None else interpolate(uexpr, V)

    return (u, uxi, uyi)

def dt_CFL(hmin, umag):
    # CFL condition: dt = 0.5 * hmin/vmax
    vmax = abs(umag.vector()[:]).max()
    return 0.5 * hmin/vmax 
    
