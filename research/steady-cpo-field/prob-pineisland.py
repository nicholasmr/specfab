#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Pine Island Galcier (PIG)
"""

from R25 import *

### Setup 

prob = R25('pineisland', p0=(-1720e3,-340e3), p1=(-1415e3,-40e3)) # make sure domain is slightly larger than meshed domain for velocity interpolation 

prob.nu_realspace_LROT = 2e-3 
prob.nu_realspace_DDRX = 5*prob.nu_realspace_LROT # this multiplier is sufficient for nonlinear solver to converge

### Run requested task

if len(sys.argv) != 2:
    raise ValueError('usage: %s [pp|solve|plot]'%(sys.argv[0]))
    
if sys.argv[1] == 'pp':
    prob.preprocess(dxy_u=2, dxy_h=2)
    prob.plot_inputs()

if sys.argv[1] == 'solve-LROT':
    FP = 'LROT'
    prob.solve(FP)
    prob.plot_solution(FP)

if sys.argv[1] == 'solve-DDRX':
    FP = 'DDRX'
    prob.solve(FP)
    prob.plot_solution(FP)

if sys.argv[1] == 'solve-altbc':
    # Test robustness of solution to changes in BCs by setting single max on inflow boundaries instead of isotropic
    FP = 'LROT'
    prob.solve(FP, altbc=True)
    prob.plot_solution(FP, suffix='-altbc')
    
if sys.argv[1] == 'plot':

    prob.kwargs_gs  = dict(left=-0.03, right=0.93, top=0.9, bottom=0.15, wspace=-0.0)
    prob.kwargs_cb  = dict(pad=0.04, aspect=22, fraction=0.04, orientation='vertical')
    prob.kwargs_lbl = dict(fontsize=11+2, frameon=False, bbox=(-0.195,1.175))
    prob.kwargs_dom = dict(columnspacing=0.8, handletextpad=0.5, handlelength=1.5, bbox_to_anchor=(1.30, 1.165))
    
    s = 1.4
    prob.figsize=(4.8*s, 2*s)

    prob.xticks_major = np.arange(-1800, -1000, 100)
    prob.xticks_minor = np.arange(-1800, -1000, 50)
    prob.yticks_major = np.arange(-400, 0, 100)
    prob.yticks_minor = np.arange(-400, 0, 50)
    prob.xlims = [-1710,-1425]
    prob.ylims = [-340,-50]

    prob.veltransform = lambda u: np.log10(u)
    prob.lvls_vel = np.arange(0.5, 3.5+0.1, 0.25)
    prob.clabel_vel = r'$\log(u)$'
    prob.lvls_strainrate = np.arange(0, 60+1e-3, 6)

    # Where to plot MODF insets
    kwargs_MODF = dict(
        axsize   = 0.155, 
        xy       = [(-1600,-230), (-1611,-275), (-1600,-295)], 
        axloc    = [(0.64+xi, 0.885) for xi in np.linspace(0, 0.205, 3)],
        lbl_bbox = (-0.5,1.4),
    )
    
    ### Make plots

    prob.plot_velocities(lbl1='a', lbl2='b')
    prob.plot_results('LROT', lbl1='c', lbl2='d', lvlsE_max=4, kwargs_MODF=kwargs_MODF)
    prob.plot_results('DDRX', lbl1='e', lbl2='f', lvlsE_max=4, kwargs_MODF=kwargs_MODF)
    prob.plot_biases('LROT')
    prob.plot_altbc_compare('LROT', lvlsE_max=4, kwargs_MODF=kwargs_MODF)
    
