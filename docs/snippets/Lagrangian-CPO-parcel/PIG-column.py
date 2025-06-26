"""
Model plug-flow CPO evolution of a Lagrangian ice column over Pine Island Glacier, Antarctica
"""

import numpy as np
from scipy.integrate import solve_ivp
import xarray

from specfabpy import specfab as sf
from specfabpy import common as sfcom
from specfabpy import plotting as sfplt

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

"""
Setup
"""

### Domain 

p0   = (-1500e3, -70e3)  # starting point (x0,y0)
xy0  = (-1710e3, -335e3) # lower left-hand  corner of region of interest
xy1  = (-1425e3, -50e3)  # upper right-hand corner of region of interest
dwin = 2 # coarsen input velocities using this moving-average window size

tend = 820 # time (in years) to trace out trajectory from staring point

### Numerics

Nt = 1000 # number of time steps taken (increase this until results are robust)
L = 12 # CPO expansion series truncation
kw_ivp = dict(method='RK45', vectorized=False) # kwargs for solve_ivp()

### CPO dynamics

# Lattice rotation
iota, zeta = 1, 0 # "deck of cards" behavior 

# DDRX
A = 1.1e7  # rate prefactor (tunable parameter)
Q = 3.36e4 # activation energy (see Richards et al. (2021) and Lilien et al. (2023))
R = 8.314  # gas constant
Gamma0 = lambda D, T: A*np.sqrt(np.einsum('ij,ji',D,D)/2)*np.exp(-Q/(R*(T+273.15))) # DDRX rate factor

"""
Determine trajectory 
"""

### Velocity field 

print('*** Loading velocity field') 

ds0 = xarray.open_mfdataset('~/ice-velocity-maps/antarctica_ice_velocity_450m_v2.nc') # MEaSUREs
uxname, uyname = 'VX', 'VY'
xr = slice(xy0[0], xy1[0])
yr = slice(xy1[1], xy0[1]) # reverse on purpose
ds = ds0.sel(x=xr,y=yr).coarsen(x=dwin, boundary='trim').mean().coarsen(y=dwin, boundary='trim').mean()

### Determine trajectory

print('*** Determining trajectory') 

def ODE_traj(t, p):
    dsi = ds.interp(x=p[0], y=p[1], method='linear')
    return [dsi[uxname].to_numpy(), dsi[uyname].to_numpy()] # [ux, uy]
        
trng = (0, tend)
ti = np.linspace(trng[0], trng[1], Nt) # time points at which to evaluate the solution
pi = solve_ivp(ODE_traj, trng, p0, t_eval=ti, vectorized=False).y # trajectory [x(t), y(t)]

### Determine velocity gradient tensor along trajectory

print('*** Determining velocity gradients') 

ds_grad = xarray.Dataset({
    'duxx': ds[uxname].differentiate('x'),
    'duxy': ds[uxname].differentiate('y'),
    'duyx': ds[uyname].differentiate('x'),
    'duyy': ds[uyname].differentiate('y'),
})
ds_gradi = ds_grad.interp(x=('points', pi[0]), y=('points', pi[1]), method='linear')

ugrad = np.zeros((Nt,3,3))
ugrad[:,0,0] = ds_gradi['duxx'].values
ugrad[:,0,1] = ds_gradi['duxy'].values
ugrad[:,1,0] = ds_gradi['duyx'].values
ugrad[:,1,1] = ds_gradi['duyy'].values
ugrad[:,2,2] = -(ugrad[:,0,0]+ugrad[:,1,1]) # incompressible ice

D = (ugrad + np.einsum('ijk->ikj',ugrad))/2 # strain rate tensor
W = (ugrad - np.einsum('ijk->ikj',ugrad))/2 # stress tensor

"""
Solve CPO evolution
"""

print('*** Solving CPO evolution')

S = D # stress tensor (assume coaxiality with strain-rate tensor; magnitude does not matter for our purpose)

#T = -30 # temperature assumed constant along flow line (deg. C)
T = -60 # no DDRX if ice is very cold

# isotropic initial state
lm, nlm_len = sf.init(L) # initialize specfab
nlm_init = np.zeros((nlm_len), dtype=np.complex64) # CPO expansion coefficients
nlm_init[0] = 1/np.sqrt(4*np.pi)

def ODE(t, nlm):
    # d/dt nlm_i = M_ij . nlm_j, where nlm_i is the state vector (aka s_i)
    I = np.argmin(np.abs(ti-t)) # index for closest point on trajectory at time t
    M  = sf.M_LROT(nlm, D[I], W[I], iota, zeta) # lattice rotation always present
    M += Gamma0(D[I],T)*sf.M_DDRX(nlm, S[I])    # DDRX active if sufficient warm
    #M += Lambda0*sf.M_CDRX(nlm)                # CDRX (neglected in this example)
    M += sf.M_REG(nlm, D[I])                    # regularization
    return np.matmul(M, nlm)
    
nlm = solve_ivp(ODE, trng, nlm_init, t_eval=ti, vectorized=False).y.T # CPO state along trajectory
mi, lami = sfcom.eigenframe(nlm) # a2 eigenvectors and eigenvalues

"""
Plot results
"""

print('*** Plotting results')

fs  = 1.0 # figure scale
fig = plt.figure(figsize=(9*fs,2.5*fs))
gs  = fig.add_gridspec(nrows=1, ncols=2, wspace=0.20)
ax1, ax2 = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])

ms = 1e-3 # mapscale (m to km)
c1,c2,c3 = 'tab:blue', 'tab:red', 'k'
ctraj = 'limegreen'

### Plot trajectory over velocity map

# Create regular grid for plotting velocity field
N = 200
xv = np.linspace(xy0[0], xy1[0], N)
yv = np.linspace(xy0[1], xy1[1], N)
X, Y = np.meshgrid(xv, yv)
dsi = ds.interp(x=xv, y=yv, method='linear') # interpolate ux and uy to new grid
UX, UY = dsi[uxname].to_numpy(), dsi[uyname].to_numpy()
U = np.sqrt(np.power(UX,2) + np.power(UY,2))

# Plot velocity map
lvls = np.logspace(0.5, 3.75, 14)
norm = colors.LogNorm(vmin=lvls[0], vmax=lvls[-1])
cs = ax1.contourf(X*ms, Y*ms, U, cmap='inferno', levels=lvls, norm=norm, extend='both')
cax = make_axes_locatable(ax1).append_axes('right', size="5%", pad=0.08)
plt.colorbar(cs, label=r'$u$ (m/yr)', cax=cax)
ax1.streamplot(X[0]*ms, Y[:,0]*ms, UX, UY, color='0.7', linewidth=0.4, density=0.8)
ax1.plot(pi[0]*ms, pi[1]*ms, '-', c=ctraj, lw=2)
ax1.plot(p0[0]*ms, p0[1]*ms, 'o', c=ctraj, markeredgewidth=2, markerfacecolor='none', markersize=10)

ax1.axis('square')
ax1.set_xlim([xy0[0]*ms, xy1[0]*ms])
ax1.set_ylim([xy0[1]*ms, xy1[1]*ms])
ax1.set_xlabel('$x$ (km)')
ax1.set_ylabel('$y$ (km)')
ax1.set_title(r'Pine Island Glacier')

### Plot eigenvalues

dl = np.sqrt(np.diff(pi[0])**2 + np.diff(pi[1])**2) # segment lengths of trajectory
l = np.concatenate([[0], np.cumsum(dl)]) # cumulative sum (distance along trajectory)

ax2.plot(l*ms, lami[:,0], c=c1, label=r'$\lambda_1$')
ax2.plot(l*ms, lami[:,1], c=c2, label=r'$\lambda_2$')
ax2.plot(l*ms, lami[:,2], c=c3, label=r'$\lambda_3$')
ax2.legend(loc='center right', fancybox=False, frameon=False)

ax2.set_ylabel(r'$\lambda_i$')
ax2.set_ylim([0,1])
ax2.set_xlabel('Distance along flow line (km)')
ax2.set_xlim([l[0]*ms,l[-1]*ms])

### Plot CPOs

# Make space for CPOs above eigenvalue plot
ax3 = make_axes_locatable(ax2).append_axes("top", size="35%", pad=0.00)
ax3.set_axis_off()

geo, prj = sfplt.getprojection(rotation=-45, inclination=50)

def plotCPO(ax, nlm, mi, p0, HW=0.23, cmap='Greys'):
    axtrans = ax.transData.transform(p0)
    trans = fig.transFigure.inverted().transform(axtrans)
    axin = plt.axes([trans[0]-HW/2, trans[1]-HW/2, HW,HW], projection=prj)
    axin.set_global()
    lvlset = [np.linspace(0.05, 0.45, 8), lambda x,p:'%.1f'%x]
    sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False, nchunk=None)
    sfplt.plotcoordaxes(axin, geo, negaxes=True, color=sfplt.c_dred, axislabels='xi')
    sfplt.plotmi(axin, mi, geo, colors=('gold',), ms=6) # eigenvectors
    return axin

# Plot CPO at selected positions along flow line
for li in np.linspace(l[-1]*0.05, l[-1]*0.9, 4):
    I = np.argmin(np.abs(l-li))
    plotCPO(ax2, nlm[I], mi[I], (li*ms,0.93))

### Save plot

plt.savefig('PIG-column.png', dpi=175, pad_inches=0.1, bbox_inches='tight')
