# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Monocrystal elastic constants (C11,C33,C55,C12,C13)
"""

import numpy as np 
import copy, sys, code # code.interact(local=locals())
from uncertainties import ufloat

### Experiment names
Cij_exprname = {\
    'Bennett1968' : 'Bennett (1968)', \
    'Gammon1983'  : 'Gammon et al. (1983)', \
    'Gagnon1988'  : 'Gagnon et al. (1988)', \
    'Dantl1968'   : 'Dantl (1968)', \
    'Brockamp1964': 'Brockamp and Querfurth (1964)', \
    'Bass1957'    : 'Bass et al. (1957)', \
    'Green1956'   : 'Green and Mackinnon (1956)', \
}

### Estimated values (see e.g. Gusmeroli et al. (2012) and Diez et al. (2015))
Cij = {\
    'Bennett1968' : [14.060e9, 15.240e9, 3.060e9, 7.150e9, 5.880e9], \
    'Gammon1983'  : [13.929e9, 15.010e9, 3.014e9, 7.082e9, 5.765e9], \
    'Green1956'   : [13.33e9, 14.28e9, 3.26e9, 6.03e9, 5.08e9], \
    'Gagnon1988'  : [14.48e9, 15.63e9, 3.15e9, 7.35e9, 5.97e9], \
    'Dantl1968'   : [13.21e9, 14.43e9, 2.89e9, 6.69e9, 5.84e9], \
    'Bass1957'    : [13.3e9, 14.2e9, 3.06e9, 6.3e9, 4.6e9], \
    'Brockamp1964': [13.63e9, 14.85e9, 3.04e9, 6.69e9, 5.15e9], \
}

### ...and their uncertainties
# percentwise uncertainties
dCijpct = {\
    'Gagnon1988': [(0.3+0.4)*1e-2, (0.2+0.4)*1e-2, (0.5+0.4)*1e-2, (0.6+0.4)*1e-2, (0.5+0.4)*1e-2], \
    'Dantl1968' : [0.3e-2, 0.4e-2, 0.7e-2, 2.0e-2, 7.0e-2], \
}
# absolute uncertainties
dCij = {\
    'Bennett1968' : [0.08e9, 0.12e9, 0.034e9, 0.15e9, 0.25e9], \
    'Gammon1983'  : [0.04e9, 0.05e9, 0.01e9, 0.04e9, 0.02e9], \
    'Green1956'   : [1.98e9, 0.54e9, 0.08e9, 0.72e9, 0.72e9], \
    'Gagnon1988'  : [dCijpct['Gagnon1988'][ii]*Cij['Gagnon1988'][ii] for ii in range(5)], \
    'Dantl1968'   : [0.04e9, 0.06e9, 0.02e9, 0.13e9, 0.41e9], \
    'Bass1957'    : [0.8e9, 0.7e9, 0.015e9, 0.8e9, 0.9e9], \
    'Brockamp1964': [0.04e9, 0.03e9, 0.015e9, 0.05e9, 0.0e9], \
}

### Estimated values as ufloats

Cij_ufloat = {\
    'Bennett1968' : [ufloat(Cij['Bennett1968'][ii],  dCij['Bennett1968'][ii])  for ii in range(5)], \
    'Gammon1983'  : [ufloat(Cij['Gammon1983'][ii],   dCij['Gammon1983'][ii])   for ii in range(5)], \
    'Green1956'   : [ufloat(Cij['Green1956'][ii],    dCij['Green1956'][ii])    for ii in range(5)], \
    'Gagnon1988'  : [ufloat(Cij['Gagnon1988'][ii],   dCij['Gagnon1988'][ii])   for ii in range(5)], \
    'Dantl1968'   : [ufloat(Cij['Dantl1968'][ii],    dCij['Dantl1968'][ii])    for ii in range(5)], \
    'Bass1957'    : [ufloat(Cij['Bass1957'][ii],     dCij['Bass1957'][ii])     for ii in range(5)], \
    'Brockamp1964': [ufloat(Cij['Brockamp1964'][ii], dCij['Brockamp1964'][ii]) for ii in range(5)], \
}

### Temperature at which the above measurements were made
Cij_T = {\
    'Bennett1968' : -10, \
    'Gammon1983'  : -16, \
    'Green1956'   : -16, \
    'Gagnon1988'  : -35.5, \
    'Dantl1968'   : -16, \
    'Bass1957'    : -16, \
    'Brockamp1964': -16, \
}

### Estimated temperature dependence

C11_Gagnon1988_T = lambda T: 136.8 - 0.289*T - 0.00178*T**2
C33_Gagnon1988_T = lambda T: 147.6 - 0.311*T - 0.00189*T**2
C55_Gagnon1988_T = lambda T:  29.7 - 0.063*T - 0.00039*T**2
C12_Gagnon1988_T = lambda T:  69.4 - 0.147*T - 0.00090*T**2
C13_Gagnon1988_T = lambda T:  56.3 - 0.119*T - 0.00073*T**2

C11_Dantl1968_T = lambda T: 12.904*(1 - 1.489e-3*T - 1.85e-6*T**2)
C33_Dantl1968_T = lambda T: 14.075*(1 - 1.629e-3*T - 2.93e-6*T**2)
C55_Dantl1968_T = lambda T:  2.819*(1 - 1.601e-3*T - 3.62e-6*T**2)
C12_Dantl1968_T = lambda T:  6.487*(1 - 2.072e-3*T - 3.62e-6*T**2)
C13_Dantl1968_T = lambda T:  5.622*(1 - 1.874e-3*T - 0.00e-6*T**2)

Cij_Tfunc = {\

    'Gammon1983': lambda T: np.array([ ufloat( (1 - 1.42e-3*(T-(-16)))*np.array(Cij['Gammon1983'])[ii], dCij['Gammon1983'][ii]) for ii in range(5)]), \
    
    'Gagnon1988': lambda T: 1e8*np.array([ \
                        ufloat(C11_Gagnon1988_T(T), dCijpct['Gagnon1988'][0]*C11_Gagnon1988_T(T)), \
                        ufloat(C33_Gagnon1988_T(T), dCijpct['Gagnon1988'][1]*C33_Gagnon1988_T(T)), \
                        ufloat(C55_Gagnon1988_T(T), dCijpct['Gagnon1988'][2]*C55_Gagnon1988_T(T)), \
                        ufloat(C12_Gagnon1988_T(T), dCijpct['Gagnon1988'][3]*C12_Gagnon1988_T(T)), \
                        ufloat(C13_Gagnon1988_T(T), dCijpct['Gagnon1988'][4]*C13_Gagnon1988_T(T)), \
                    ]), \
    
    'Dantl1968':  lambda T: 1e9*np.array([ \
                        ufloat(C11_Dantl1968_T(T), dCijpct['Dantl1968'][0]*C11_Dantl1968_T(T)), \
                        ufloat(C33_Dantl1968_T(T), dCijpct['Dantl1968'][1]*C33_Dantl1968_T(T)), \
                        ufloat(C55_Dantl1968_T(T), dCijpct['Dantl1968'][2]*C55_Dantl1968_T(T)), \
                        ufloat(C12_Dantl1968_T(T), dCijpct['Dantl1968'][3]*C12_Dantl1968_T(T)), \
                        ufloat(C13_Dantl1968_T(T), dCijpct['Dantl1968'][4]*C13_Dantl1968_T(T)), \
                    ]), \
}

### Inferred values

g_exprname = {\
    '003_0': r'\#003, $\alpha=0$',\
    '003_1': r'\#003, $\alpha=1$',\
    '003_5': r'\#003, $\alpha=1/2$',\
    '007_0': r'\#007, $\alpha=0$',\
    '007_1': r'\#007, $\alpha=1$',\
    '007_5': r'\#007, $\alpha=1/2$',\
    '010_0': r'\#010, $\alpha=0$',\
    '010_1': r'\#010, $\alpha=1$',\
    '010_5': r'\#010, $\alpha=1/2$',\
}

g_inferred = {\
	'003_0': [6.0726e+00, 3.3188e+00, 0.8684, 0.8106, 1.1170],\
	'003_1': [5.9903e+00, 3.3137e+00, 0.9084, 0.7929, 1.0949],\
	'003_5': [6.0234e+00, 3.3177e+00, 0.8913, 0.8010, 1.1045],\
	'007_0': [6.4994e+00, 3.3860e+00, 0.8335, 0.8713, 1.0612],\
	'007_1': [6.4334e+00, 3.3872e+00, 0.8557, 0.8610, 1.0552],\
	'007_5': [6.4638e+00, 3.3870e+00, 0.8453, 0.8658, 1.0580],\
	'010_0': [6.3397e+00, 3.3192e+00, 0.8431, 0.8265, 1.1072],\
	'010_1': [6.2391e+00, 3.3203e+00, 0.8831, 0.8107, 1.0898],\
	'010_5': [6.2825e+00, 3.3206e+00, 0.8655, 0.8180, 1.0974],\
}

### Functions

def Cij_to_g(Cij):
    C11,C33,C55,C12,C13 = Cij # unpack
    C66 = (C11-C12)/2
    lam  = C12
    mu   = C66
    gam  = lam + 2*mu
    Elam = C13/C12
    Emu  = C55/C66
    Egam = C33/C11
    g = [lam,mu,gam,Elam,Emu,Egam] # note: contains \gamma which is not independent 
    return np.array(g)

