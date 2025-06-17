#!/usr/bin/python3
# N. Rathmann <rathmann@nbi.ku.dk>, 2024-

import numpy as np
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from steadyCPOfancy import steadyCPOfancy

### Domain and plotting

domain = dict(
    name        = 'ross',
    fmeasures   = '~/ice-velocity-maps/antarctica_ice_velocity_450m_v2.nc',
    fbedmachine = '~/ice-velocity-maps/BedMachineAntarctica-v3.nc',
    fgeo        = 'mesh.geo',
    subsample_u = 2,
    subsample_h = 4,
)

scpo = steadyCPOfancy(domain)

# Setup plotting options
scpo.figsize = (8.2,3)
scpo.kw_gs   = dict(wspace=0.25)
scpo.kw_cax  = dict(size="3.5%", pad=0.11)
scpo.kw_leg['bbox_to_anchor'] = (-0.04, 1.2)
scpo.xticks_major = np.arange(-1100, 600, 400)
scpo.xticks_minor = np.arange(-1100, 600, 200)
scpo.yticks_major = np.arange(-1400, 0, 300)
scpo.yticks_minor = np.arange(-1400, 0, 150)
scpo.kw_vel['kw_tcf']['levels'] = np.logspace(0.5, 3.0, 11)
scpo.kw_epsE['kw_tcf']['levels'] = np.arange(0, 6+1e-3, 1)
scpo.kw_E['kw_tcf']['levels'] = np.arange(0, 3+.01, 0.25)
scpo.kw_E['kw_cb']['ticks'] = scpo.kw_E['kw_tcf']['levels'][::4]
scpo.kw_MODF = dict(
        axsize   = 0.145, 
        xy       = [(-180,-1080),(-20,-1050),(180,-1100),(-70,-680)], 
        axloc    = [(0.6+xi, 0.76) for xi in np.linspace(0, 0.235, 4)],
        lbl_bbox = (-0.55,1.4),
)

### Problem specification

problem_LROT  = dict(name='LROT',  T=None, bcs=[[1,scpo.bc_isotropic],])
problem_altbc = dict(name='altbc', T=None, bcs=[[1,scpo.bc_zsinglemax],])
problem_DDRX  = dict(name='LROT+DDRX', T=np.linspace(-40, -15, 4), bcs=problem_LROT['bcs'])
numerics      = dict(L=8, nu_orimul=0.8, nu_real=2e-3, nu_realmul=[50,20,12,8])

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

