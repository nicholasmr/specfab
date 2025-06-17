#!/usr/bin/python3
# N. Rathmann <rathmann@nbi.ku.dk>, 2024-

import numpy as np
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from steadyCPOfancy import steadyCPOfancy

### Domain and plotting

domain = dict(
    name        = 'pineisland',
    fmeasures   = '~/ice-velocity-maps/antarctica_ice_velocity_450m_v2.nc',
    fbedmachine = '~/ice-velocity-maps/BedMachineAntarctica-v3.nc',
    fgeo        = 'mesh.geo',
    subsample_u = 2,
    subsample_h = 2,
)

scpo = steadyCPOfancy(domain)

# Setup plotting options
scpo.figsize = (7.2,3)
scpo.kw_gs   = dict(wspace=0.41)
scpo.kw_cax  = dict(size="4.5%", pad=0.11)
scpo.kw_leg['bbox_to_anchor'] = (-0.07, 1.17)
scpo.xticks_major = np.arange(-2000, -1000, 100)
scpo.xticks_minor = np.arange(-2000, -1000, 50)
scpo.yticks_major = np.arange(-500, 200, 100)
scpo.yticks_minor = np.arange(-500, 200, 50)
scpo.kw_vel['kw_tcf']['levels'] = np.logspace(0.5, 3.5, 13)
scpo.kw_epsE['kw_tcf']['levels'] = np.arange(0, 50+.01, 5)
scpo.kw_MODF = dict(
    axsize   = 0.145, 
    xy       = [(-1600,-230), (-1611,-275), (-1600,-295)], 
    axloc    = [(0.63+xi, 0.835) for xi in np.linspace(0, 0.19, 3)],
    lbl_bbox = (-0.5,1.4),
)

### Problem specification

problem_LROT  = dict(name='LROT', T=None, bcs=[[1,scpo.bc_isotropic],])
problem_altbc = dict(name='altbc', T=None, bcs=[[1,scpo.bc_zsinglemax],])
problem_DDRX  = dict(name = 'LROT+DDRX', T=np.linspace(-40, -15, 4), bcs=problem_LROT['bcs'])
numerics      = dict(L=8, nu_orimul=0.8, nu_real=3e-4, nu_realmul=[50,20,12,8])

### Run requested task

if len(sys.argv) != 2:
    print('usage: %s [pp|solve-LROT|solve-DDRX|solve-altbc]'%(sys.argv[0]))
    sys.exit(1)
    
if sys.argv[1] == 'pp':
    scpo.preprocess()
    scpo.plot_inputs()
            
if sys.argv[1] == 'solve-LROT':
    scpo.solve(problem_LROT, numerics)
    scpo.plot_results(problem_LROT, numerics)

if sys.argv[1] == 'solve-DDRX':
    scpo.solve(problem_DDRX, numerics)
    scpo.plot_results(problem_DDRX, numerics)

if sys.argv[1] == 'solve-altbc':
    scpo.solve(problem_altbc, numerics)
    scpo.plot_results(problem_altbc, numerics)
    
if sys.argv[1] == 'plot-biases':
    problem = problem_LROT
    scpo.femsolution(problem['name'], numerics['L'])
    scpo.plot_biases(problem)

