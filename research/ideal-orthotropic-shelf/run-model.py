#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

import sys, os, code # code.interact(local=locals())
import numpy as np
from datetime import datetime

from lib.plottools import *
from lib.model import *
from R25config import R25config

### Initialize flow model

if len(sys.argv) != 3:
    print('usage: %s [LROT|DDRX] L'%(sys.argv[0]))
    sys.exit(0)

conf = R25config() 
conf.rheology = "Orthotropic"
#conf.rheology = "Isotropic"  # debug
conf.fabdyn   = str(sys.argv[1]) # fabric evolution type
conf.L        = int(sys.argv[2]) # spectral space truncation

model = FlowModel(conf.Nt, **conf.kwargs_flowmodel())

os.system('mkdir -p %s'%(conf.path_diagnostics()))

### Integrate model

tstart = datetime.now()
tt_save = 0

for tt in np.arange(conf.nt, conf.Nt):

    dt = 1e-10 if tt == 0 else model.dt_CFL() # time step size
    
    print("\n==================================")
    print('Step %i of %i :: dt = %.3e'%(tt, conf.Nt, dt))
    print("==================================")
    print("Run time: %s\n"%(str(datetime.now()-tstart)))
    
    # Modify Newton solver parameters for Stokes problem
#    tolmul = 1e-5 if tt==conf.nt else 1e+1 # VERY SAFE
    tolmul = 1e-5 if tt==conf.nt else 5e+2 # GOOD PERFORMANCE WITH THIS SETTING
#    tolmul = 5e3 # DEBUG (fast first step)

    relaxation = 0.40 if tt==conf.nt else 0.40
    
    model.integrate(dt, tolmul=tolmul, relaxation=relaxation)
    
    if tt % conf.tt_savemod == 0:
        plot_diagnostics(model, conf, conf.fname_diagn(tt_save))
        model.save(conf.fname_state(tt_save))
        tt_save += 1
    
