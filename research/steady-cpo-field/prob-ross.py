#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Ross
"""

from R25 import *

### Setup 

prob = R25('ross', p0=(-920e3,-1420e3), p1=(620e3,-380e3)) # make sure domain is slightly larger than meshed domain for velocity interpolation

prob.nu_realspace_LROT = 2e-3
prob.nu_realspace_DDRX = 3*prob.nu_realspace_LROT # this multiplier is sufficient for nonlinear solver to converge

### Run requested task

if len(sys.argv) != 2:
    raise ValueError('usage: %s [pp|solve|plot]'%(sys.argv[0]))
    
if sys.argv[1] == 'pp':
    prob.preprocess(dxy_u=2, dxy_h=8)
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

    prob.kwargs_gs  = dict(left=0.095, right=0.93, top=1.03, bottom=0.03, wspace=0.25, hspace=0.08)
    prob.kwargs_cb  = dict(pad=0.05, aspect=18.5, fraction=0.032, orientation='vertical')
    prob.kwargs_lbl = dict(fontsize=11+2, frameon=False, bbox=(-0.175,1.21))
    prob.kwargs_dom = dict(columnspacing=0.8, handletextpad=0.5, handlelength=1.5, bbox_to_anchor=(1.05, 1.2))

    s = 1.1
    prob.figsize=(7*s, 1/3*6.8*s)

    prob.xticks_major = np.arange(prob.x0km, prob.x1km+1, 300)
    prob.xticks_minor = np.arange(prob.x0km, prob.x1km+1, 150)
    prob.yticks_major = np.arange(-1500, -200+1, 200)
    prob.yticks_minor = np.arange(-1500, -200+1, 100)
    prob.xlims = [-900, 600]
    prob.ylims = [-1335, -400]

    prob.veltransform = lambda u: np.log10(u)
    prob.lvls_vel = np.arange(0.3,3.1+1e-3, 0.2)
    prob.clabel_vel = r'$\log(u)$'
    prob.lvls_strainrate = np.arange(0, 6+1e-3, 1)

    # Where to plot MODF insets
    kwargs_MODF = dict(
        axsize   = 0.18, 
        xy       = [(-180,-1080),(-20,-1050),(180,-1100),(-70,-680)], 
        axloc    = [(0.625+xi, 0.86) for xi in np.linspace(0, 0.25, 4)],
        lbl_bbox = (-0.50,1.35),
    )

    ### Make plots

    prob.plot_velocities(lbl1='a', lbl2='b')
    prob.plot_results('LROT', lbl1='c', lbl2='d', kwargs_MODF=kwargs_MODF)
    prob.plot_results('DDRX', lbl1='e', lbl2='f', kwargs_MODF=kwargs_MODF)
    prob.plot_biases('LROT')
    prob.plot_altbc_compare('LROT', kwargs_MODF=kwargs_MODF)

