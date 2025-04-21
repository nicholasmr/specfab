# N. M. Rathmann <rathmann@nbi.ku.dk>, 2024

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import time

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, rc

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt
FS = 11
FS = sfplt.setfont_tex(fontsize=FS)
FSAX = FS+0

#---------------------
# Setup
#---------------------

### Experiments

L_list = (8,10,12)
NN = 15
Ngrains_list = np.int64(np.logspace(2, np.log10(3000), NN))

t_SDM = np.zeros(len(L_list))
t_DDM = np.zeros(len(Ngrains_list))

### Deformation kinematics

Nparcels = 13500
Nt = 2000 

axis  = 2 # axis of shortening (T>0) or lengthening (T<0): 0=x, 1=y, 2=z
Tc    = 1 # time taken in seconds for parcel to reduce to half (50%) height if T>0, or abs(time) taken for parcel to double in height (200%) if T<0.
r     = 0 # asymmetry parameter for shortening (if T>0) or lengthening (if T<0)

T = Tc/np.log(2) # corresponding e-folding time (tau)
ugrad = sf.pureshear_ugrad(axis, r, T) # velocity gradient
D, W = sf.ugrad_to_D_and_W(ugrad)      # strain-rate and spin tensor

t = 4 # some specific time of interest
b   = sf.pureshear_r(T, t)          # scaling parameter at time t
F   = sf.pureshear_F(axis, r, T, t) # deformation gradient tensor at time t
eps = sf.F_to_strain(F)             # strain tensor at time t

print(eps)

dt = t/Nt
Dn = np.repeat(D, Nt, axis=0)
Wn = np.repeat(W, Nt, axis=0)

#---------------------
# Integrate 
#---------------------

for ii,L in enumerate(L_list):
    lm, nlm_len = sf.init(L)
    nlm0 = np.zeros((nlm_len), dtype=np.complex64)
    t_start = time.time()
    nlm = sf.nlm_LROT(nlm0, dt, Nt, Dn, Wn, +1)
    t_SDM[ii] = time.time()-t_start

lm, nlm_len = sf.init(8) # L value does not matter here
for ii,N in enumerate(Ngrains_list):
    r = np.array([1,0,0], dtype=np.float32) # does not matter all axes are the same
    ri0 = np.tile(r.T, (N,1))
    print(ri0.shape)
    t_start = time.time()
    ri = sf.ri_LROT(ri0, dt, Nt, Dn, Wn, +1)
    t_DDM[ii] = time.time()-t_start

# We only model one of the slip-system axes here, so if two are modelled (as required for olivine), the time would double
t_SDM *= 2 * Nparcels/Nt
t_DDM *= 2 * Nparcels/Nt

Ndofs_list = [int((L+2)**2/4) for L in L_list]
size_SDM = np.array([Nparcels *       2 * sys.getsizeof(np.array([1e0]*Ndofs)) for Ndofs   in Ndofs_list])
size_DDM = np.array([Nparcels * Ngrains * sys.getsizeof(np.array([1e0]*3))     for Ngrains in Ngrains_list]) # three euler angles per grain per parcel

#---------------------
# Plot
#---------------------

scale=0.8
fig = plt.figure(figsize=(8*scale,2.7*scale))
gs = gridspec.GridSpec(1, 2, wspace=0.35, hspace=0, bottom=0.22, top=0.86, left=0.1, right=0.98)

cdred = '#a50f15'
ls = ['-','--',':']

### Speed

ax = fig.add_subplot(gs[0,0])

ax.plot(Ngrains_list, t_DDM, '-', c='k', label='DDM')
for ii,L in enumerate(L_list):
    ax.plot([0,5e3], [t_SDM[ii]]*2, ls[ii], c=cdred, label='$L=%i$'%(L))
    
ax.set_xlabel('Number of grains')
ax.set_ylabel("Integration time for\n %i parcels/nodes (s)"%(Nparcels))
kwargs_leg = dict(fancybox=False, frameon=False, ncols=2, columnspacing=1, handlelength=1.2, labelspacing=.2)
ax.legend(loc=2, bbox_to_anchor=(0,1.03), **kwargs_leg)

ax.set_xscale('log')
ax.set_xlim(Ngrains_list[[0,-1]])
ax.set_yticks(np.arange(0,70,10))
ax.set_yticks(np.arange(0,70,5), minor=True)
ax.set_ylim([0,40])

sfplt.panellabel(ax, 2, r'\textit{(a)}', frameon=False, fontsize=FS+2, bbox=(-0.165,1.3))

### Memory

ax = fig.add_subplot(gs[0,1])

scale = 1e-6 # mb
ax.plot(Ngrains_list, size_DDM*scale, '-', c='k', label='DDM/D-Rex')
for ii,L in enumerate(L_list):
    ax.plot([0,5e3], [size_SDM[ii]*scale]*2, ls[ii], c=cdred, label='$L=%i$'%(L))
    
ax.set_xlabel('Number of grains')
ax.set_ylabel("Memory required for\n %i parcels/nodes (MB)"%(Nparcels))
ax.legend(loc=5, bbox_to_anchor=(1,0.35), **kwargs_leg)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(Ngrains_list[[0,-1]])

sfplt.panellabel(ax, 2, r'\textit{(b)}', frameon=False, fontsize=FS+2, bbox=(-0.185,1.3))

### FINISH

fout = 'SDM-DDM-performance.pdf'
print('Saving %s'%(fout))
#fig.set_tight_layout(True) 
plt.savefig(fout, dpi=200)
plt.close()
            
