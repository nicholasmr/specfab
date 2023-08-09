# The f2py fortran module specfabpy.so provides "specfabpy", which is the interface we actually want to provide the user with.
# To avoid this double (nested) naming, we import the fortran interface directly to the main namespace.
from .specfabpy import specfabpy as specfab 

