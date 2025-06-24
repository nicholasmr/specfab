"""
Modeled CPO profile of the GRIP ice core, Greenland
"""

import numpy as np
from scipy import interpolate
import pandas as pd

import matplotlib.pyplot as plt

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from specfabpy import specfab as sf
from specfabpy import common as sfcom
from specfabpy import plotting as sfplt

### Init

L = 12 # expansion series truncation
lm, nlm_len = sf.init(L)

### Velocity gradient experienced by parcel 

H = 3027 # ice thickness (Montagnat et al., 2014)
a = 0.24 # meter ice equiv. per yr (Montagnat et al., 2014)

tau = H/a # e-folding time scale
ugrad = -1/tau * np.diag([-0.5, -0.5, 1]) # uniaxial compression along z-axis
D = (ugrad+np.transpose(ugrad))/2 # symmetric part (strain rate tensor)
W = (ugrad-np.transpose(ugrad))/2 # anti-symmetric part (spin tensor)
S = D # stress tensor (assume coaxiality with strain-rate tensor; magnitude does not matter for our purpose)

### Fabric dynamics

# Lattice rotation
iota, zeta = 1, 0 # "deck of cards" behavior 

# DDRX
A = 1.1e7 # rate prefactor (tunable parameter)
Q = 3.36e4 # activation energy (see Richards et al. (2021) and Lilien et al. (2023))
R = 8.314  # gas constant
Gamma0 = lambda D, T: A*np.sqrt(np.einsum('ij,ji',D,D)/2)*np.exp(-Q/(R*(T+273.15))) # DDRX rate factor

### Numerics 

Nt = 500 # number of time steps
dt = 100 # time step size (yr)
ti = np.arange(0,Nt) * dt # time vector
zi = np.exp(-ti/tau) # relative height above bed at each point in time

### Temperature profile 

df = pd.read_csv('../../../data/icecores/GRIP/temperature.csv') # fetch from github
f = interpolate.interp1d(df['zrel'].to_numpy(), df['T'].to_numpy(), kind='nearest', fill_value='extrapolate')
Ti = f(zi) # temperature vector
#Ti[:] = -60 # no DDRX if very cold
    
### Initial fabric state

nlm  = np.zeros((Nt,nlm_len), dtype=np.complex64) # state vector
lami = np.zeros((Nt,3)) # a2 eigenvalues

lxy = 0.25 # initial horizontal eigenvalues
a2_0 = np.diag([lxy, lxy, 1-2*lxy]) # initial a2 surface state
nlm[0,:sf.L2len] = sf.a2_to_nlm(a2_0) # initial state vector
lami[0] = sfcom.eigenframe(nlm[0])[1] # eigenvalues of initial state

### Euler integration

for tt in np.arange(1,Nt):
    nlm_0 = nlm[tt-1,:] # previous solution
    T = Ti[tt] # temperature from borehole measurements
    M_LROT = sf.M_LROT(nlm_0, D, W, iota, zeta) # lattice rotation operator
    M_DDRX = Gamma0(D,T)*sf.M_DDRX(nlm_0, S)    # DDRX operator
    M_REG  = sf.M_REG(nlm_0, D)                 # regularization operator
    M      = M_LROT + M_DDRX + M_REG
    nlm[tt]  = nlm_0 + dt*np.matmul(M, nlm_0) # Euler step
    lami[tt] = sfcom.eigenframe(nlm[tt])[1]

### Plot modeled eigenvalues

fig = plt.figure(figsize=(3,4))
ax = plt.subplot(111)

c1,c2,c3 = 'tab:green', 'tab:red', 'k'

ax.plot(lami[:,0], zi, '-',  c=c1, label=r'$\lambda_1$')
ax.plot(lami[:,1], zi, '-',  c=c2, label=r'$\lambda_2$')
ax.plot(lami[:,2], zi, '--', c=c3, label=r'$\lambda_3$')

ax.legend(loc=1, fancybox=False, frameon=False)
ax.set_title(r'GRIP ice core')

ax.set_xlabel(r'$\lambda_i$')
ax.set_xticks(np.arange(0,1+.01,0.2))
ax.set_xlim([0,1])

ax.set_ylabel(r'$z/H$')
ax.set_yticks(np.arange(0,1+.01,0.1))
ax.set_ylim([0,1])

### Plot modeled CPOs

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
        
for _ in np.linspace(0.1, 0.9, 4):
    I = np.argmin(np.abs(zi-_))
    plotCPO(ax, nlm[I], (1.2,_))

### Plot observations 

df = pd.read_csv('../../../data/icecores/GRIP/orientations.csv') # fetch from github
zi = df['zrel'].to_numpy()

kw = dict(marker='o', facecolor='none', zorder=1)
ax.scatter(df['lam1'].to_numpy(), zi, edgecolor=c1, **kw)
ax.scatter(df['lam2'].to_numpy(), zi, edgecolor=c2, **kw)
ax.scatter(df['lam3'].to_numpy(), zi, edgecolor=c3, **kw)

### Save plot

plt.savefig('GRIP.png', dpi=175, pad_inches=0.1, bbox_inches='tight')
