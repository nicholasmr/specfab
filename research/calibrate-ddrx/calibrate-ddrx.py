# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien

"""
Calibrate DDRX activation function energy (Q) and prefactor (A) by reproducing Dome C (EDC) fabric profile.
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from scipy.interpolate import interp1d
import pickle
from progress.bar import Bar
import matplotlib.pyplot as plt

from specfabpy import specfab as sf
from specfabpy import common as sfcom
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

yr2s = 31556926
s2yr = 3.16887646e-8
c2k  = lambda degc: degc+273.15 # deg. C to deg. K

#---------------------
# Numerics
#---------------------

L = 12   # spectral truncation
Nt = 100 # number of integration steps

#---------------------
# Gamma_0(T) params
#---------------------

# Guess 
Q = 0.5e5 # activation energy
A = 2e10  # prefactor

print('Gamma(T=-10)/Gamma(T=-20) = %e'%(sf.Gamma0(np.eye(3),c2k(-10),A,Q)/sf.Gamma0(np.eye(3),c2k(-20),A,Q)))

#---------------------
# Load EDC profile
#---------------------

file = open("EDC.p",'rb')
(fab,temp,aux) = pickle.load(file)
n20_0 = 0.14 # initial fabric state of parcel (only l,m=2,0 component is assumed nonzero)
strain_zz_stop=-0.96 # stop integration at this vertical parcel strain. Fabric is only measured until the depth corresponding to this vertical strain.

### Geometry

z_bed = -aux['H'] # depth of bed
z0 = abs(fab['z'][0]) 
z1 = aux['H']
H = z1-z0 # ice column height as seen by parcel model
print('z_bed, H = %f, %f'%(z_bed, H))

### Measured fabric

lami_meas, z_meas = fab['eig'], fab['z']

### Temperature interpolant

T, zT = temp['T'], temp['z']
f_T = interp1d(zT,T) # C. Ritz (pers. comm.)

### Accumulation rate

#b = aux['b'] # present day
b = 0.0153 # average over core (D. Lilien)
print('b = %3e m/yr'%(b))

#---------------------
# Ice parcel model
#---------------------
    
### Mode of deformation
# Assumes depth-constant vertical uniaxial compression: the classical Nye model of a dome.

b  /= yr2s
t_e = 1/(b/H) # e-folding time for uniaxial compression (strainrate_zz = MeanAccum/H)
r   = 0 # uniaxial compression
ax  = 2 # z axis
    
### Constants

tend    = sf.pureshear_strainii_to_t(strain_zz_stop, t_e)
timevec = np.linspace(0,tend,Nt)
dt      = timevec[1]-timevec[0] # this is the step size needed if we want Nt time steps given t_e
strain_zz = np.array([sf.F_to_strain(sf.pureshear_F(ax,r,t_e,t))[-1,-1] for ii, t in enumerate(timevec)]) # for the constant time-step size, these are the vertical parcel strains as a function of time
z = H*(strain_zz+1) - H # Vertical parcel strains correspond to these depths
z += -z0 # effective offset of parcel model

### Initialize problem 

lm, nlm_len = sf.init(L)

nlm = np.zeros((Nt, nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[0,0] = 1/np.sqrt(4*np.pi) # normalized distribution
nlm[0,3] = n20_0 # initial single max strength (manually fitted to match eigenvalue profile)

lami = np.zeros((Nt, 3)) # array of modelled eigenvalues
lami[0,:] = np.diag(sf.a2(nlm[0,:])) # a2 already diagonal for this mode of deformation, so simply extract the diagonal values

D, W = sf.ugrad_to_D_and_W(sf.pureshear_ugrad(ax,r,t_e)) # strain-rate and spin tensors for mode of deformation
T = f_T(z) # temperature at modelled parcel depths 

gam = np.zeros((Nt)) # DDRX rate factor magnitude
gam[0] = sf.Gamma0(D, c2k(T[0]), A, Q)

### Euler integration of fabric model

with Bar('dt=%.3fyr, Nt=%i :: L=%i (nlm_len=%i) ::'%(dt*s2yr,Nt,L,nlm_len), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds') as bar:
    for nn in np.arange(1,Nt):
    
        nlm_prev = nlm[nn-1,:]
        
        ### Lattice rotation
        
        M_LROT = sf.M_LROT(nlm_prev, D, W, 1, 0) 

        ### DDRX
        
        S = D.copy() # assume stress and strain-rate are coaxial (not necessarily true)
        gam[nn] = sf.Gamma0(S, c2k(T[nn]), A, Q) # DDRX decay rate magnitude
        M_DDRX = gam[nn]*sf.M_DDRX(nlm_prev, S) # DDRX operator

        ### Regularization
        
        M_REG = sf.M_REG(nlm_prev, D) 

        ### All together
        
        M = M_LROT + M_DDRX + M_REG
        nlm[nn,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # Euler step
        lami[nn,:] = np.diag(sf.a2(nlm[nn,:]))
        
        bar.next() # update progress bar
        
# Estimate numerical accuracy 
diff_abs = np.abs(100*(nlm[:,0]-nlm[0,0])/nlm[0,0])
print(r'Numerical error estimate: max(n00(t) - n00(t=0)) [%] = ', np.amax(diff_abs))

#--------------
# Plot results
#--------------

scale = 0.55
fig = plt.figure(figsize=(14*scale,10*scale))
gs  = fig.add_gridspec(ncols=4, nrows=3, width_ratios=[1.8,1.3,1,1], right=0.97, top=0.95, hspace=0.55)
ax1 = fig.add_subplot(gs[:, 0])
ax2 = fig.add_subplot(gs[:, 2], sharey=ax1)
ax3 = fig.add_subplot(gs[:, 3], sharey=ax1)

### Eigenvalues

ax1.scatter(lami_meas[:,0], z_meas, c=sfplt.c_lgreen)
ax1.scatter(lami_meas[:,1], z_meas, c=sfplt.c_lblue)
ax1.scatter(lami_meas[:,2], z_meas, c=sfplt.c_lred)

ax1.plot(lami[:,0],z, c=sfplt.c_dgreen, ls='-',  lw=2, clip_on=False)
ax1.plot(lami[:,1],z, c=sfplt.c_dblue,  ls='--', lw=2, clip_on=False)
ax1.plot(lami[:,2],z, c=sfplt.c_dred,   ls='-',  lw=2, clip_on=False)

ax1.set_xlabel(r'$\lambda_i$')
ax1.set_xlim([0,1])
ax1.set_ylabel(r'$z$ (m)')
ax1.set_ylim([z_bed,0])
ax1.set_yticks(-np.flipud(np.arange(0,H,500)))
ax1.set_yticks(-np.flipud(np.arange(0,H,250)), minor=True)

ax1.set_title(r'{\bf Dome C} // $A=%.1e$, $Q=%.1e$ // $L=%i$'%(A,Q,L), fontsize=FS)

### Temperature

#ax2.plot(f_T(zT), zT, 'g', lw=2)
ax2.plot(T, z, 'k', lw=2)
ax2.set_xlabel(r'$T$ ($\SI{}{\degreeCelsius}$)')
plt.setp(ax2.get_yticklabels(), visible=False)

### Gamma

ax3.semilogx(gam, z, 'k', lw=2)
ax3.set_xlabel(r'$\Gamma_0$')
plt.setp(ax3.get_yticklabels(), visible=False)
ax3.set_xticks([1e-16,1e-14,1e-12])

### ODFs

geo, prj = sfplt.getprojection(rotation=45, inclination=45)
axi = [fig.add_subplot(gs[_, 1], projection=prj) for _ in range(3)]
for _ in axi: _.set_global() 

lvlset = (np.linspace(0.0, 0.7, 8), lambda x,p:'%.1f'%x)
for ii, nn in enumerate([int(Nt*2/10), int(Nt*2/3), -1]):
    ax1.plot([0,1],[z[nn],]*2, ':', c='k', lw=2)
    sfplt.plotODF(nlm[nn,:], lm, axi[ii], lvlset=lvlset, cblabel='ODF at $z$=%im'%(z[nn]))
    sfplt.plotcoordaxes(axi[ii], geo, axislabels='vuxi')
    
fname = 'calibrate-ddrx-L%i-%s.png'%(L, 'EDC')
print('saving %s'%(fname))
plt.savefig(fname, dpi=150)

