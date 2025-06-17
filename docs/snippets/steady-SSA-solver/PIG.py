"""
Example of steady state SSA fabric (CPO) solver for Pine Island Glacier, Antarctica
"""

import numpy as np
from specfabpy.fenics.steadyCPO import steadyCPO

"""
Model domain
"""

domain = dict(
    name        = 'PIG', # All output file names be prepended by this domain name
    fmeasures   = '~/ice-velocity-maps/antarctica_ice_velocity_450m_v2.nc', # MEaSUREs ice velocities
    fbedmachine = '~/ice-velocity-maps/BedMachineAntarctica-v3.nc',         # BedMachine
    fgeo        = 'mesh.geo', # Gmsh finite element mesh
    subsample_u = 2, # Coarsen input velocities using a moving average with window size "subsample_u" (i.e., no subsampling unless > 1)
    subsample_h = 2, # ...same but for ice thicknesses
)

scpo = steadyCPO(domain)

### Generate mesh and interpolate ice velocities and ice geometry onto mesh
scpo.preprocess() # Once run, this can be outcommented

"""
Steady CPO problem
"""

problem = dict(
    ### Solve for: lattice rotation + advection
    # - If lattice rotation is dominant then set the temperature to T=None (negligible DDRX; cold ice limit)
    name = 'LROT', 
    T    = None, 
    
    ### Solve for: lattice rotation + DDRX + advection
    # - If DDRX cannot be neglected, the rate of DDRX is calculated provided the ice temperature and strain rate fields
    # - The rate of DDRX is parameterized using the lab-calibrated temperature activation function of Lilien et al. (2023) (can be changed)
    # - The temperature field T is assumed constant (can be changed)
    # - The solver works by gradually approaching the nonlinear LROT+DDRX solution by starting out assuming no DDRX, then 
    #   gradually increasing the temperature in discrete steps until reaching the desired temperature 
#    name = 'LROT+DDRX', 
#    T    = np.linspace(-40, -15, 4), # Steps of increasing temperatures used to converge on solution for -15C
    
    ### Boundary conditions (BCs)
    # - BCs are given as a list of mesh "Physical Line" IDs (see mesh.geo) and the CPO state to be imposed there
    # - Two different CPO states can be specified on boundaries: isotropic (scpo.bc_isotropic) or a perfect vertical single maximum (scpo.bc_zsinglemax)
    # - If no BC is specified on a given boundary, it is assumed free (unconstrained)
    bcs = [[1,scpo.bc_isotropic],] # Isotropic ice on boundary ID=1, unconstrained elsewhere
)

numerics = dict(
    ### Spectral truncation
    # - For most cases L>=8 should be sufficient, but setting L>=12 might require too many resources on laptop hardware
    L = 8,
    
    ### S^2 (orientation space) regularization multiplier
    # - Multiplier of the strength of the default orientation space (S^2) regularization
    # - In general, this *should* be nu_orimul >= 0.8
    nu_orimul  = 0.8,
    
    ### R^3 (real space) regularization strength
    # - Free parameter that *must* be adjusted for the specific problem/domain
    # - If set too small the solver will hang (not converge), whereas if too large the CPO field will become unphysically smooth
    nu_real = 0.3e-3, 

    ### R^3 (real space) regularization multiplier *for DDRX-activated problems*
    # - List of successively weaker R^3 regularization multipliers, used to converge on the DDRX-activated solution
    # - DDRX requires in general stronger regularization (i.e., nu_realmul > 1)
    nu_realmul = [50,20,12,8] 
)

### Run the solver and dump solution (may take a few minutes depending on problem size)
scpo.solve(problem, numerics) # Once run, this can be outcommented

"""
Plot results
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

ms2myr = 3.17098e+8
m2km = 1e-3
mapscale = m2km # x and y axis scale

def newfig(boundaries=True, floating=True, mesh=False, bgcolor='0.85', figsize=(3,3)):
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    legh, legt = [], []
    if bgcolor is not None: 
        x0, x1 = mapscale*(scpo.x0-scpo.dx), mapscale*(scpo.x1+scpo.dx)
        y0, y1 = mapscale*(scpo.y0-scpo.dy), mapscale*(scpo.y1+scpo.dy)
        ax.add_patch(plt.Rectangle((x0,y0), x1-x0, y1-y0, color=bgcolor))
    if mesh: 
        ax.triplot(triang, lw=0.075, color='0.5', alpha=0.8, zorder=10)
    if boundaries:
        coords, *_ = scpo.bmesh(problem['bcs'])
        colors = ['cyan', 'magenta', 'yellow']
        markers = ['s',]*len(colors)
        for ii, (xb, yb) in enumerate(coords):
            ax.scatter(xb*mapscale, yb*mapscale, c=colors[ii], marker=markers[ii], s=3, zorder=12, clip_on=False)
            legh.append(Line2D([0], [0], color=colors[ii], lw=2))
            legt.append('Isotropic' if ii == 0 else 'Free')
    if floating: 
        ax.tricontour(triang, mask==3, [0.5, 1.5], colors=['limegreen',], linewidths=2, zorder=11)
        legh.append(Line2D([0], [0], color='limegreen', lw=2))
        legt.append('Floating')
        
    ax.legend(legh, legt, ncol=3, loc=1, bbox_to_anchor=(1.15,1.15), fancybox=False, frameon=False, \
                handletextpad=0.5, columnspacing=0.8, handlelength=1.3)
    ax.axis('square')
    ax.set_xlabel(r'$x$ (km)')
    ax.set_ylabel(r'$y$ (km)')
    ax.set_xlim([mapscale*scpo.x0, mapscale*scpo.x1])
    ax.set_ylim([mapscale*(scpo.y0-scpo.dy), mapscale*scpo.y1])
    return (fig, ax)
    
def newcax(ax): return make_axes_locatable(ax).append_axes("right", size="4%", pad=0.13)

kw_save = dict(dpi=150, pad_inches=0.1, bbox_inches='tight')

### Load saved fields 

coords,cells, ux,uy,umag,epsE, S,B,H,mask = scpo.npinputs()
coords,cells, mi,lami,E_CAFFE = scpo.npsolution(problem['name'])
triang = scpo.triang(coords, cells, mapscale=mapscale)

### Velocities

fig, ax = newfig(boundaries=False)
lvls = np.logspace(0.5, 3.5, 13)
cs = ax.tricontourf(triang, ms2myr*umag, levels=lvls, norm=colors.LogNorm(vmin=lvls[0], vmax=lvls[-1]), extend='both', cmap='inferno')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label('$u$ (m/yr)')
plt.savefig('%s-umag.png'%(domain['name']), **kw_save)

### Effective strain rate

fig, ax = newfig(boundaries=False)
lvls = np.arange(0, 50+.01, 5)
cs = ax.tricontourf(triang, 1e3*ms2myr*epsE, levels=lvls, extend='max', cmap='viridis')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$\dot{\epsilon}_{e}$ (1/yr)')
plt.savefig('%s-epsE.png'%(domain['name']), **kw_save)

### Delta lambda (horizontal eigenvalue difference)

fig, ax = newfig()
lvls = np.arange(0, 0.8+.01, 0.1)
dlam = abs(lami[:,0] - lami[:,1]) # eigenvalues 1 and 2 are the largest and smallest in-model-plane eigenvalues
cs = ax.tricontourf(triang, dlam, levels=lvls, extend='max', cmap='Blues')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$\Delta\lambda$')

# Quiver principal horizontal eigenvector
meshpts = (coords[0,:], coords[1,:])
xv, yv = np.linspace(scpo.x0, scpo.x1, 15)[1:-1], np.linspace(scpo.y0, scpo.y1, 15)[1:-1]
x, y = np.meshgrid(xv, yv, indexing='xy')
m1 = mi[:,0,:] # principal horizontal eigenvector
m1x = griddata(meshpts, m1[:,0].flatten(), (x, y), method='linear', fill_value=np.nan)
m1y = griddata(meshpts, m1[:,2].flatten(), (x, y), method='linear', fill_value=np.nan) # y coordinate is index 2 (z coordinate) since problem is in xz plane
renorm = np.sqrt(m1x**2+m1y**2)
m1x, m1y = np.divide(m1x, renorm), np.divide(m1y, renorm)
hq = ax.quiver(mapscale*x, mapscale*y, +m1x, +m1y, color='tab:red', scale=40)
hq = ax.quiver(mapscale*x, mapscale*y, -m1x, -m1y, color='tab:red', scale=40)
ax.quiverkey(hq, 0.1, 0.05, 3, r'${\bf m}_1$', labelpos='E')

plt.savefig('%s-%s-dlam.png'%(domain['name'], problem['name']), **kw_save)

### lambda_z (vertical eigenvalue)

fig, ax = newfig()
lvls = np.arange(0, 0.8+.01, 0.1)
lamz = lami[:,2] # eigenvalue 3 is the out-of-model-plane (z) eigenvalue
cs = ax.tricontourf(triang, lamz, levels=lvls, extend='max', cmap='RdPu')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$\lambda_z$')
plt.savefig('%s-%s-lamz.png'%(domain['name'], problem['name']), **kw_save)

### E (CAFFE)

fig, ax = newfig()
lvls = np.logspace(-1, 1, 17)
cs = ax.tricontourf(triang, E_CAFFE, levels=lvls, norm=colors.LogNorm(vmin=lvls[0], vmax=lvls[-1]), extend='both', cmap='PuOr_r')
hcb = plt.colorbar(cs, cax=newcax(ax))
hcb.set_label(r'$E$', labelpad=-2)
plt.savefig('%s-%s-E.png'%(domain['name'], problem['name']), **kw_save)

