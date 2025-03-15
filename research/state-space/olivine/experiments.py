# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

"""
Experiment definitions.

Samples should be approximately rotationally symmetric for this analysis to be meaningful.
If not, correlations may not follow modeled line, or incorrectly seem to fall outside the range of physically-allowed eigenvalues.
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import glob

expr_Boneh2014 = {\
    'path': 'Boneh_etal_2014', \
    'plotname': 'Boneh and Skemer (2014)', \
    'type':'uc', \
    'source':'experiment', \
    'coords':'spherical', 'I_azi':1, 'I_lat':0, 'isdeg':False, 'iscolat':True,
#    'flist': ['116.csv.csv', '114.csv.csv','125.csv.csv','122.csv.csv','126.csv.csv', ], \
#    'Fzz': np.exp(-np.array([0.217,0.328,0.406,0.507,0.653])), 
#    'T': [1200]*5, \
    'flist': ['114.csv','125.csv','122.csv','126.csv', ], \
    'Fzz': np.exp(-np.array([0.328,0.406,0.507,0.653])), 
    'T': [1200]*4, \
    'skiprows': 0, \
    'sep': ',', \
    'color':'0.4' \
}

expr_Kumamoto2019 = {\
    'path': 'Kumamoto_etal_2019', \
    'plotname': 'Kumamoto et al. (2019)', \
    'type':'ss', \
    'source':'natural', \
    'coords':'spherical', 'I_azi':1, 'I_lat':0, 'isdeg':False, 'iscolat':True,
    'flist': ['3923J08.ctf.csv', '3923J09.ctf.csv', '3923J11.ctf.csv', '3924J03a.ctf.csv', '3924J08.ctf.csv', '3924J09a.ctf.csv', 'JP10M11.ctf.csv'], \
    'Fzz': 1+np.array([0.81, 0.81, 0.65, 3.37, 5.25, 3.86, 1.36]), \
    'T': [0]*7, \
    'skiprows': 1, \
    'sep': ',', \
    'color':'tab:green' \
}

expr_Bernard2019 = {\
    'path': 'Bernard_etal_2019', \
    'plotname': 'Bernard et al. (2019)', \
    'type':'ss', \
    'source':'natural', \
    'coords':'spherical', 'I_azi':1, 'I_lat':0, 'isdeg':False, 'iscolat':True,
    'flist': ['07EB4-01.crc.csv', '114027-6.crc.csv', '114027-16.crc.csv', '114027-23.crc.csv', 'BELB9-6a.crc.csv', 'CG07-1-52.crc.csv', \
              'DW8.crc.csv', 'N122.crc.csv', 'O18A.crc.csv', 'R-DH-26.crc.csv', 'WCiVb47.crc.csv' ], \
    'skiprows': 1, \
    'sep': ',', \
    'color':'tab:pink' \
}

expr_Yabe2020 = {\
    'path': 'Yabe_etal_2020', \
    'plotname': 'Yabe et al. (2020)', \
    'type':'uc', \
    'source':'experiment', \
    'coords':'spherical', 'I_azi':1, 'I_lat':0, 'isdeg':False, 'iscolat':True,
    'flist': ['KH71_221207_OrientationData.ang.csv',], \
    'Fzz': [0.49,], # h0/h = exp(-0.7) \
    'T': [1250,], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:purple' \
}

expr_Kim2022 = {\
    'path': 'Kim_etal_2022', \
    'plotname': 'Kim et al. (2022)', \
    'type':'cc', \
    'source':'experiment', \
    'coords':'spherical', 'I_azi':1, 'I_lat':0, 'isdeg':False, 'iscolat':True,
    'flist': ['R-KH%i.ctf.csv'%(ii) for ii in [186,201,211,214,215,216,231,233,236,237,238]], \
    'skiprows': 1, \
    'sep': ',', \
    'color':'tab:blue' \
}


expr_Miyazaki2013_uc = {\
    'path': 'Miyazaki_etal_2013__UC', \
    'plotname': 'Miyazaki et al. (2013)', \
    'type':'uc', \
    'source':'experiment', \
#    'coords':'vector', \
    'coords':'spherical', 'I_azi':1, 'I_lat':0, 'isdeg':False, 'iscolat':True,
    'flist': ['FoDi20_com1250_KF-157.ctf.csv', 'FoDi20_com1350_KF-191.ctf.csv'], \
    'Fzz': [0.6, 0.6], # or should these be =exp(-0.6) like Yabe et al.? \
    'T': [1250, 1350], # deg. C \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:cyan' \
}

expr_Miyazaki2013_ue = {\
    'path': 'Miyazaki_etal_2013__UE', \
    'plotname': 'Miyazaki et al. (2013)', \
    'type':'ue', \
    'source':'experiment', \
#    'coords':'vector', \
    'coords':'spherical', 'I_azi':1, 'I_lat':0, 'isdeg':False, 'iscolat':True,
    'flist': ['FoDi20_ten1250_KS-10.ctf.csv', 'FoDi20_ten1350_KS-14.ctf.csv'], \
    'Fzz': [1+0.6, 1+0.6], \
    'T': [1250, 1350], # deg. C \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:orange' \
}

#--------------
# D-REX model trajectories
#--------------

dn = 10

I = np.arange(0, 2*250 +1, dn)

expr_DREX_uc = {\
    'path': '../drex/state-trajectories/axisymmetricCompression', \
    'plotname': 'drex UC', \
    'type':'uc', \
    'source':'drex', \
    'coords':'vector', \
    'flist': ['drex-tstep%i.csv'%(ii) for ii in I], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:red' \
}

expr_DREX_ue = {\
    'path': '../drex/state-trajectories/axisymmetricExtension', \
    'plotname': 'drex UE', \
    'type':'ue', \
    'source':'drex', \
    'coords':'vector', \
    'flist': ['drex-tstep%i.csv'%(ii) for ii in I], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:red' \
}

I = np.arange(0, 3*250 +1, dn)

expr_DREX_ss = {\
    'path': '../drex/state-trajectories/simpleShear', \
    'plotname': 'drex SS', \
    'type':'ss', \
    'source':'drex', \
    'coords':'vector', \
    'flist': ['drex-tstep%i.csv'%(ii) for ii in I], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:red' \
}

