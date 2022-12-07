# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Experiment definitions.

Samples should be approximately rotationally symmetric for this analysis to be meaningful.
If not, correlations may not follow modeled line, or incorrectly seem to fall outside the range of physically-allowed eigenvalues.
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import glob

expr_Hunter = {\
    'path': 'Hunter_etal', \
    'plotname': 'H22, unconf. compr., warm', \
    'type':'ucw', \
    'coords':'quaternion', \
    'flist': ['MD14/md14.csv','MD12/md12.csv','MD4/md4.csv','MD22/md22.csv', \
              'MD9/md9.csv','MD3/md3.csv','MD13/md13.csv', \
              'D1-1/D1-1_pf.csv','D1-5/D1-5_pf.csv','D1-7/D1-7_pf.csv', \
              'D5-1/D5-1_pf.csv','D5-3/D5-3_pf.csv'], \
#    'flist': ['D5-1/D5-1_pf.csv','D5-3/D5-3_pf.csv'], # D5-1 and D5-3 were added above for Lilien's work \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:gray' \
}

# (Qi et al., 2019)
expr_Qi = {\
    'path': 'Qi_etal', \
    'plotname': 'Q19, simple shear, cold', \
    'type':'ss', \
    'coords':'quaternion', \
    'flist': ['PIL135 shear 20 clean.ctf.csv', 'PIL142 shear 20 clean.ctf.csv', 'PIL144_shear_30_clean.ctf.csv', 'PIL145_shear_30_clean.ctf.csv'], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:olive' \
}
            
# (Fan et al., 2020, 2021)
expr_Fan_30 = {\
    'path': 'Fan_etal_-30', \
    'plotname': 'F20, unconf. compr., cold', \
    'type':'uc', \
    'coords':'quaternion', \
    'flist': ['series_-30/PIL165_30mu.ctf.csv', 'series_-30/PIL166_30mu.ctf.csv', 'series_-30/PIL243_30mu.ctf.csv'], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'#6a3d9a' \
}

expr_Fan_20 = {\
    'path': 'Fan_etal_-20', \
    'plotname': None, \
    'type':'uc', \
    'coords':'quaternion', \
    'flist': ['series_-20/PIL182_30mu.ctf.csv', 'series_-20/PIL184_30mu.ctf.csv', 'series_-20/PIL185_30mu.ctf.csv', 'series_-20/PIL254_30mu.ctf.csv', 'series_-20/PIL255_30mu.ctf.csv'], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'#6a3d9a' \
}

expr_Fan_10 = {\
    'path': 'Fan_etal_-10', \
    'plotname': 'F20, unconf. compr., warm', \
    'type':'ucw', \
    'coords':'quaternion', \
#    'flist': ['series_-10/PIL007_30mu.crc.csv','series_-10/PIL163_30mu.crc.csv','series_-10/PIL176_30mu.crc.csv','series_-10/PIL177_30mu.crc.csv','series_-10/PIL178_30mu.crc.csv',], \
    'flist': ['series_-10/PIL007_30mu.crc.csv','series_-10/PIL163_30mu.crc.csv','series_-10/PIL177_30mu.crc.csv','series_-10/PIL178_30mu.crc.csv',], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'#e31a1c' \
}

expr_Fan_4 = {\
    'path': 'Fan_etal_-4', \
    'plotname': None, \
    'type':'ucw', \
    'coords':'quaternion', \
#    'flist': ['series_-4/undeformed.ctf.csv','series_-4/OIL006.ctf.csv','series_-4/OIL007.ctf.csv','series_-4/OIL008.ctf.csv','series_-4/OIL009.ctf.csv',], \
    'flist': ['series_-4/OIL006.ctf.csv','series_-4/OIL007.ctf.csv','series_-4/OIL008.ctf.csv'], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'#e31a1c' \
}

expr_Priestley = {\
    'path': 'Priestley', \
    'plotname': 'Priestley shear margin', \
    'type':'ss', \
    'coords':'quaternion', \
    'flist': ['003.ctf.csv', '007.ctf.csv', '013.ctf.csv', '023.ctf.csv', '036.ctf.csv', '047.ctf.csv', '063.ctf.csv'], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'#33a02c' \
}

files = np.sort([int(os.path.basename(f)[0:4]) for f in glob.glob("data/SPICE/*SPICE/*.txt")]) #.sort(key=int) #.sort(key=natural_keys)
expr_SPICE = {\
    'path': 'SPICE', \
    'plotname': 'SPICE', \
    'type':'ue', \
    'coords':'spherical', 'I_azi':2, 'I_lat':1, 'iscolat':True, \
    'flist': ["%i SPICE/%i SPICE.txt"%(fid,fid) for fid in files[::1]], \
    'skiprows': 22, \
    'sep': '\t', \
    'color':'#ff7f00' \
}

# UNPUBLISHED DATA; Contact Ilka Weikusat (AWI) for permission to use
expr_EGRIP_S5 = {\
    'path': 'EGRIP_S5', \
    'plotname': 'EGRIP shear margin', \
    'type':'ss', \
    'coords':'spherical', 'I_azi':0, 'I_lat':1, 'iscolat':False, \
    'flist': ['stereo_EGRIP_S5_124_1_20.txt', 'stereo_EGRIP_S5_124_3_20.txt', 'stereo_EGRIP_S5_124_volume_vertical_2_20.txt'], \
    'skiprows': 0, \
    'sep': '\t', \
    'color':'#6a3d9a' \
}

# UNPUBLISHED DATA; Contact Ilka Weikusat (AWI) for permission to use
expr_EGRIP_MAIN = {\
    'path': 'EGRIP_MAIN', \
    'plotname': 'EGRIP', \
    'type':'ue', \
    'coords':'spherical', 'I_azi':0, 'I_lat':1, 'iscolat':False, \
    'flist': ['EGRIP3328_1_stereo.txt', 'EGRIP3680_3_stereo.txt', 'EGRIP3795_3_stereo.txt'], \
    'skiprows': 0, \
    'sep': '\t', \
    'color':'#6a3d9a' \
}

expr_GRIP = {\
    'path': 'GRIP', \
    'plotname': 'GRIP', \
    'type':'uc', \
    'coords':'spherical', 'I_azi':0, 'I_lat':1, 'iscolat':False, \
    'flist': ['SPLITTED/GRIP_CAXES_%02i'%(ii) for ii in np.arange(1,28)], \
    'skiprows': 3, \
    'sep': '\t', \
    'color':'#1f78b4' \
}

files = np.sort([int(os.path.basename(f)[3:]) for f in glob.glob("LAWDOME/SPLITTED/DSS*")]) #.sort(key=int) #.sort(key=natural_keys)
expr_LAWDOME = {\
    'path': 'LAWDOME', \
    'plotname': 'Law Dome', \
    'type':'uc', \
    'coords':'spherical', 'I_azi':-1, 'I_lat':-2, 'iscolat':True, \
#    'flist': ["SPLITTED/DSS%i"%(fid) for fid in files[::1]], \
    'flist': ["SPLITTED/DSS%i"%(fid) for fid in [116, 182, 199, 214, 254, 266, 287, 326, 356, 391, 465, 473, 481, 496, 521, 542, 604, 662, 697, 732, 796, 852, 888, 902, 930, 971, 1005, 1017, 1049, 1099]], \
    'skiprows': 0, \
    'sep': ',', \
    'color':'tab:pink' \
}

