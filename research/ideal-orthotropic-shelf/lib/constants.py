#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024

"""
Physical constants, conversion factors, etc.
"""

### Physical constants
GRAV_ACCEL    = 9.8 # m/s
DENSITY_ICE   = 917 # kg/m^3
DENSITY_OCEAN = 1024 # kg/m^3
GAS_CONST     = 8.314462618 # J/(K*mol)

### Conversion factors
MS_TO_MYR  = 31556926
MS_TO_KMYR = MS_TO_MYR/1000
MYR_TO_MS  = 1/MS_TO_MYR

### Boundary/domain IDs used by fenics
DOM_ID__IN  = 1
DOM_ID__TOP = 2
DOM_ID__OUT = 3
DOM_ID__BOT = 4

