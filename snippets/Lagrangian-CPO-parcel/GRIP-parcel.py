"""
Model CPO profile of the GRIP ice core, Greenland
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy import interpolate
import pandas as pd

from specfabpy import specfab as sf
from specfabpy import common as sfcom
from specfabpy import plotting as sfplt

import matplotlib.pyplot as plt

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

"""
Setup
"""

### Velocity gradient experienced by parcel 

H   = 3027 # ice thickness (Montagnat et al., 2014)
a   = 0.24 # meter ice equiv. per yr (Montagnat et al., 2014)
tau = H/a  # e-folding time scale
ugrad = -1/tau*np.diag([-0.5, -0.5, 1]) # uniaxial compression along z-axis

### Numerics 

tend = 50e3 # time (in years) to trace out trajectory, starting from the surface
#tend = -tau*np.log(0.05) # alternatively, set tend so that simulation stops at 95% thinning

Nt = 1000 # number of time steps taken (increase this until results are robust)
ti = np.linspace(0, tend, Nt) # time points at which to evaluate the solution
z  = np.exp(-ti/tau) # relative height above bed at each point in time    

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

### Temperature profile 

df = pd.read_csv('../../../data/icecores/GRIP/temperature.csv') # fetch from github
fz = interpolate.interp1d(df['zrel'].to_numpy(), df['T'].to_numpy(), kind='linear', fill_value='extrapolate')
T = fz(z) # temperature profile
    
""" 
Solve CPO evolution
"""

D = (ugrad+np.transpose(ugrad))/2 # symmetric part (strain rate tensor)
W = (ugrad-np.transpose(ugrad))/2 # anti-symmetric part (spin tensor)
S = D # stress tensor (assume coaxiality with strain-rate tensor; magnitude does not matter for our purpose)

#T[:] = -60 # no DDRX if ice is very cold

lm, nlm_len = sf.init(L) # initialize specfab    
lh_init = 0.25 # initial horizontal eigenvalues
a2_init = np.diag([lh_init, lh_init, 1-2*lh_init]) # initial a2 (at surface)

nlm_init = np.zeros((nlm_len), dtype=np.complex64) # CPO expansion coefficients
nlm_init[:sf.L2len] = sf.a2_to_nlm(a2_init) # initial state vector (at surface)

def ODE(t, nlm):
    # d/dt nlm_i = M_ij . nlm_j, where nlm_i is the state vector (aka s_i)
    I = np.argmin(np.abs(ti-t)) # index for closest point on trajectory at time t
    M  = sf.M_LROT(nlm, D, W, iota, zeta) # lattice rotation always present
    M += Gamma0(D,T[I])*sf.M_DDRX(nlm, S) # DDRX active if sufficient warm
    #M += Lambda0*sf.M_CDRX(nlm)          # CDRX (neglected in this example)
    M += sf.M_REG(nlm, D)                 # regularization
    return np.matmul(M, nlm)

nlm = solve_ivp(ODE, (0, tend), nlm_init, t_eval=ti, vectorized=False).y.T # CPO state along trajectory
mi, lami = sfcom.eigenframe(nlm) # a2 eigenvectors and eigenvalues

"""
Plot results
"""

### Plot modeled eigenvalues

fig = plt.figure(figsize=(3,4))
ax = plt.subplot(111)

c1,c2,c3 = 'tab:blue', 'tab:red', 'k'

ax.plot(lami[:,0], z, '-',  c=c1, label=r'$\lambda_1$')
ax.plot(lami[:,1], z, '-',  c=c2, label=r'$\lambda_2$')
ax.plot(lami[:,2], z, '--', c=c3, label=r'$\lambda_3$')

ax.legend(loc=1, fancybox=False, frameon=False)
ax.set_title(r'GRIP ice core')

ax.set_xlabel(r'$\lambda_i$')
ax.set_xticks(np.arange(0,1+.01,0.2))
ax.set_xlim([0,1])

ax.set_ylabel(r'$z/H$')
ax.set_yticks(np.arange(0,1+.01,0.1))
ax.set_ylim([0,1])

### Plot CPOs

geo, prj = sfplt.getprojection(rotation=45, inclination=50)

def plotCPO(ax, nlm, p0, HW=0.2, cmap='Greys'):
    axtrans = ax.transData.transform(p0)
    trans = fig.transFigure.inverted().transform(axtrans)
    axin = plt.axes([trans[0]-HW/2, trans[1]-HW/2, HW,HW], projection=prj)
    axin.set_global()
    lvlset = [np.linspace(0.05, 0.45, 8), lambda x,p:'%.1f'%x]
    sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False, nchunk=None)
    sfplt.plotcoordaxes(axin, geo, negaxes=False, color=sfplt.c_dred, axislabels='xi')
    return axin
        
for zi in np.linspace(0.1, 0.9, 4):
    I = np.argmin(np.abs(z-zi))
    plotCPO(ax, nlm[I], (1.2,zi))

### Plot observations 

df = pd.read_csv('../../../data/icecores/GRIP/orientations.csv') # fetch from github
zobs = df['zrel'].to_numpy()
kw = dict(marker='o', facecolor='none', zorder=1)
ax.scatter(df['lam1'].to_numpy(), zobs, edgecolor=c1, **kw)
ax.scatter(df['lam2'].to_numpy(), zobs, edgecolor=c2, **kw)
ax.scatter(df['lam3'].to_numpy(), zobs, edgecolor=c3, **kw)

### Save plot

plt.savefig('GRIP-parcel.png', dpi=175, pad_inches=0.1, bbox_inches='tight')
