# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

"""
Fabric state-space diagram comparing model trajectories to data (Lilien et al., 2022)
"""

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import quaternion as qt # pip3 install numpy-quatern
import pandas as pd
import pickle, glob
from progress.bar import Bar

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)
sys.path.insert(0, '../../../demo') # for importing local specfabpy build (if available) and common python header
from header import * # contains matplotlib setup etc.
from specfabpy import specfabpy as sf

SELFNAME = sys.argv[0][:-3] # used as prefix for pickled files
os.system('mkdir -p specfab-state-trajectories')
def pfile(fname): return "specfab-state-trajectories/%s--%s.p"%(SELFNAME, fname) # full path name for pickled files

#--------------------
# Flags
#--------------------

INTEGRATE_MODEL = 1  # Generate model lines? Else load saved Pickle files.
DEBUG           = 0  # Low-resolution plotting, etc.
            
#--------------------
# Config
#--------------------

Nt = 300 # Number of specfab integration steps

if DEBUG: 
    L = 8 # Spectral truncation for debugging (faster plotting)
    RES = 2*50 
    RESX, RESY = int(RES/2), RES
else:
    L = 20 # Spectral truncation 
    RESX = RESY = 5*100
    
#--------------------
# Experimental data to plot
#--------------------

# All data
#experiments = (expr_GRIP, expr_LAWDOME,  expr_EGRIP_MAIN, expr_SPICE, expr_EGRIP_S5, expr_Priestley, expr_Qi, expr_Fan_10,expr_Fan_4, expr_Hunter) 

# Lilien et al. (2022)
experiments = (expr_GRIP, expr_LAWDOME,  expr_EGRIP_MAIN, expr_SPICE, expr_Priestley, expr_Qi, expr_Fan_10,expr_Fan_4, expr_Hunter) 

# DEBUG
#experiments = (expr_GRIP,)

#--------------------
# Load pre-processed correlations (see process-experiments.py)
#--------------------

correlations = []

for expr in experiments:
    fcorr = "observed-states/%s.p"%(expr['path']) 
    print('=== Loading correlations: %s ==='%(fcorr))
    corr = pickle.load(open(fcorr, "rb"))
    corr['n20'] /= normfac # Normalization
    corr['n40'] /= normfac
    corr['n60'] /= normfac
    correlations.append(corr) 
    
#--------------------
# Modeled correlations
#--------------------

lm, nlm_len = sf.init(L) 

def integrate_model(nlm0, Mtype, ugrad, dt, Nt=Nt, rotate=False, name=None):

    eps = (ugrad+np.transpose(ugrad))/2 # Symmetric part (strain-rate)
    omg = (ugrad-np.transpose(ugrad))/2 # Anti-symmetric part (spin)

    nlm = np.zeros((Nt,nlm_len), dtype=np.complex64)
    nlm[0,:] = nlm0 # Initial state 
    nlmr = nlm.copy() # nlm array for rotated c-axis distribution
             
    # Euler integration 
    with Bar('dt=%.3e, Nt=%i :: L=%i (nlm_len=%i) ::'%(dt,Nt,L,nlm_len), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds') as bar:
        for tt in np.arange(1,Nt):

            nlm_prev = nlm[tt-1,:]
                
            if Mtype == 'LROT':
                iota, zeta = 1, 0
                M_LROT = sf.M_LROT(nlm_prev, eps, omg, iota, zeta)
                M_REG  = sf.M_REG(nlm_prev, eps)
                M      = M_LROT + M_REG
                
            if Mtype == 'DDRX':
                M = 1e-0*sf.M_DDRX(nlm_prev,eps)

            if Mtype == 'CDRX':
                M = 1e-2*sf.M_CDRX(nlm_prev)

            nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev)

            if rotate:
                a2 = sf.a2(nlm[tt,:])
                (v1_colat, v1_lon, _) = get_v_angles(a2)
                nlmr[tt,:] = sf.rotate_nlm( nlm[tt,:], 0, -v1_lon)
                nlmr[tt,:] = sf.rotate_nlm(nlmr[tt,:], +v1_colat, 0)
            else:
                nlmr[tt,:] = nlm[tt,:] # no rotation requested, just copy the unrotated solution
            
            bar.next()

    if name is not None:
        pickle.dump([nlm, nlmr, lm, nlm_len], open(pfile(name), "wb"))

    return nlm, nlmr


if INTEGRATE_MODEL:

    print('*** Generating model trajectories from scratch. You can re-use the saves trajectories to avoid re-calculating them by setting INTEGRATE_MODEL=1')

    ### Solve for state-vector time evolution
    
    nlm_iso = np.zeros((nlm_len), dtype=np.complex64)
    nlm_iso[0] = normfac
    nlm_cc, _  = integrate_model(nlm_iso, 'LROT', +np.diag([.5, .5, -1]), 0.014, name='cc') # confined compression (cc)
    nlm_ce, _  = integrate_model(nlm_iso, 'LROT', -np.diag([.5, .5, -1]), 0.007, name='ce') # confined extension (ce)
    _, nlmr_ss = integrate_model(nlm_iso, 'LROT', +np.array([[0,1,0], [0,0,0], [0,0,0]]), 0.030, name='ss', rotate=True) # simple vertical shear (ss)

    f_nlm_init = lambda nlm, n20_ref: nlm[np.argmin(np.abs(nlm[:,3]-n20_ref)), :] 
    nlm_ddrx1, _ = integrate_model(f_nlm_init(nlm_ce, -0.20), 'DDRX', +np.diag([.5, .5, -1]), 5*0.05, name='ddrx1') # DDRX trajectory 1
    nlm_ddrx2, _ = integrate_model(f_nlm_init(nlm_cc, +0.00), 'DDRX', +np.diag([.5, .5, -1]), 5*0.05, name='ddrx2') # DDRX trajectory 2
    nlm_ddrx3, _ = integrate_model(f_nlm_init(nlm_cc, +0.20), 'DDRX', +np.diag([.5, .5, -1]), 5*0.05, name='ddrx3') # DDRX trajectory 3
    nlm_ddrx4, _ = integrate_model(f_nlm_init(nlm_cc, +0.40), 'DDRX', +np.diag([.5, .5, -1]), 5*0.05, name='ddrx4') # DDRX trajectory 4
    
    nlm_cdrx1, _ = integrate_model(f_nlm_init(nlm_ce, -0.27), 'CDRX', +np.diag([.5, .5, -1]), 5*0.05, name='cdrx1') # DDRX trajectory 1
    nlm_cdrx2, _ = integrate_model(f_nlm_init(nlm_cc, +0.32), 'CDRX', +np.diag([.5, .5, -1]), 5*0.05, name='cdrx2') # DDRX trajectory 1
    
### Load solutions
nlm_cc, _,  lm, nlm_len = pickle.load(open(pfile('cc'), "rb"))
nlm_ce, _,  lm, nlm_len = pickle.load(open(pfile('ce'), "rb"))
_, nlmr_ss, lm, nlm_len = pickle.load(open(pfile('ss'), "rb"))

nlm_ddrx1, _, lm, nlm_len = pickle.load(open(pfile('ddrx1'), "rb"))
nlm_ddrx2, _, lm, nlm_len = pickle.load(open(pfile('ddrx2'), "rb"))
nlm_ddrx3, _, lm, nlm_len = pickle.load(open(pfile('ddrx3'), "rb"))    
nlm_ddrx4, _, lm, nlm_len = pickle.load(open(pfile('ddrx4'), "rb"))

nlm_cdrx1, _, lm, nlm_len = pickle.load(open(pfile('cdrx1'), "rb"))
nlm_cdrx2, _, lm, nlm_len = pickle.load(open(pfile('cdrx2'), "rb"))
                
### Normalize
nlm_cc    = np.array([ nlm_cc[tt,:]/nlm_cc[tt,0]   for tt in np.arange(Nt) ])
nlm_ce    = np.array([ nlm_ce[tt,:]/nlm_ce[tt,0]   for tt in np.arange(Nt) ])
nlmr_ss   = np.array([ nlmr_ss[tt,:]/nlmr_ss[tt,0] for tt in np.arange(Nt) ])
nlm_ddrx1 = np.array([ nlm_ddrx1[tt,:]/nlm_ddrx1[tt,0] for tt in np.arange(Nt) ])
nlm_ddrx2 = np.array([ nlm_ddrx2[tt,:]/nlm_ddrx2[tt,0] for tt in np.arange(Nt) ])
nlm_ddrx3 = np.array([ nlm_ddrx3[tt,:]/nlm_ddrx3[tt,0] for tt in np.arange(Nt) ])
nlm_ddrx4 = np.array([ nlm_ddrx4[tt,:]/nlm_ddrx4[tt,0] for tt in np.arange(Nt) ])
nlm_cdrx1 = np.array([ nlm_cdrx1[tt,:]/nlm_cdrx1[tt,0] for tt in np.arange(Nt) ])
nlm_cdrx2 = np.array([ nlm_cdrx2[tt,:]/nlm_cdrx2[tt,0] for tt in np.arange(Nt) ])


#--------------------
# Polynomial-fitted correlation curve used by Gerber et al. (2022)
#--------------------

if 0:
    ### ***NEEDS UPDATING IF TO BE USED AGAIN***
    
    Ne_girdle = 0
    Ne_smax   = 1
    x_model = np.hstack(( np.flipud(np.real(nlmr_sf[:,Ne_girdle,3])),  np.real(nlmr_sf[:,Ne_smax,3])))
    y_model = np.hstack(( np.flipud(np.real(nlmr_sf[:,Ne_girdle,10])), np.real(nlmr_sf[:,Ne_smax,10])))

    pc = np.polyfit(x_model, y_model, 4)  # Last argument is degree of polynomial
    #print("Coefficient values:\n", pc)
    p_model = np.poly1d(pc)
    p_fit = p_model(x_model)
    #print(p_fit-y_model)

    pickle.dump([x_model, y_model, pc, p_model], open("corrpoly.p", "wb"))


#--------------------
# Construct plot
#--------------------

ms = 6.0
mse = 7 # end-member case points
FSLEG = FS-1.5
FSANNO = FS-2.5

c_girdle = '#8c510a' 
c_smax   = '#01665e'
c_ddrx   = '#b2182b'
c_cdrx   = 'k'

cl_smax   = np.array([199,234,229])/255
cl_girdle = np.array([246,232,195])/255
cl_circle = np.array([253,219,199])/255

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

### Setup figure

scale = 4.7
fig = plt.figure(figsize=(1.8*scale,1.1*scale))
plt.subplots_adjust(left=0.08, right=0.725, top=0.90, bottom=0.125)
ax = plt.gca()
xlims, ylims = [-1.45,2.65], [-1.8,3.75]
sc = np.diff(ylims)/np.diff(xlims)
x = np.linspace(xlims[0],xlims[1],RESX)
y = np.linspace(ylims[0],ylims[1],RESY)

### Determine valid subspace (valid eigenvalues)

validregion = np.zeros((RESY, RESX)) # 0 = invalid, 1 = valid 
validregion_lowerbound = np.zeros((RESX)) # points along lower boundary, used for colouring the background (shading) of subspace with ~circle fabrics.
print('Determining subspace of valid eigenvalues...', end='')
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        nlm_ = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients
#        nlm_[0], nlm_[3], nlm_[10] = normfac, x_, y_
        nlm_[0], nlm_[3], nlm_[10] = 1, x_, y_
        a2_ = sf.a2(nlm_) # diagional
        a2_eigvals = np.sort(np.diag(a2_))
        Q1,Q2,Q3,Q4,Q5,Q6, a4_eigvals = sf.a4_eigentensors(nlm_)
        isvalid_a2_eig = (np.amin(a2_eigvals)>=0) and (np.amax(a2_eigvals)<=1)  
        isvalid_a4_eig = (np.amin(a4_eigvals)>=0) and (np.amax(a4_eigvals)<=1)
        validregion[yii,xii] = isvalid_a2_eig and isvalid_a4_eig
    validregion_lowerbound[xii] = y[np.argmax(validregion[:,xii])]
            
print('done')

### Determine subspace shadings

imdat = np.empty((RESY, RESX, 4), dtype=float) # color (0,1,2) and alpha (3)

# LATROT
xmodel_latrot = np.concatenate((nlm_ce[:, 3],nlm_cc[:, 3]))
ymodel_latrot = np.concatenate((nlm_ce[:,10],nlm_cc[:,10]))
I_latrot = np.argwhere(np.abs(xmodel_latrot)>0.1/normfac) # make white shading near isotropy
xmodel_latrot = xmodel_latrot[I_latrot]
ymodel_latrot = ymodel_latrot[I_latrot] 

# DDRX
I_DDRX = np.arange(int(RESY*3.3/10),int(RESX*5.8/10))
xmodel_DDRX = x[I_DDRX]
ymodel_DDRX = validregion_lowerbound[I_DDRX]

# Join all model points
xmodel = np.append(xmodel_latrot, xmodel_DDRX)
ymodel = np.append(ymodel_latrot, ymodel_DDRX)

print('Determining shading regions for different fabric types...', end='')
y_cutoff = np.interp(x, (xlims[0],0.1/normfac,xlims[-1]), np.array([-0.15,-0.15,0.125])/normfac)
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        if validregion[yii,xii]:
            distnn = np.amin( np.sqrt(np.real((x_-xmodel)**2 + (1/sc*(y_-ymodel))**2)) ) # distance from model-line points
            var, expo, ycutoff = 1e-3, 6, -0.14/normfac
            imdat[yii,xii,-1] = np.exp(-distnn**expo/var) # set alpha depending on distance to model lines/points
            if y_ > y_cutoff[xii]: # assume single max or girdle fabric
                if x_<0: imdat[yii,xii,0:-1] = cl_girdle
                else:    imdat[yii,xii,0:-1] = cl_smax
            else: # circle fabric
               imdat[yii,xii,0:-1] = cl_circle
               imdat[yii,xii,-1] = np.exp(-distnn**expo/(1.75*var))
        else: 
            imdat[yii,xii,0:-1] = [0.85]*3 # bad 
            imdat[yii,xii,-1] = 1
print('done')

# Plot valid/invalid subspaces and fabric-type shadings within valid subspace
im = ax.imshow(imdat, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)

### Determine E_zz

Ezz = np.zeros((RESY, RESX)) 
Ecc,Eca,alpha,nprime = 1, 1e3, 0.0125, 1
print('Determining E_zz ...', end='')
for xii, x_ in enumerate(x):
    for yii, y_ in enumerate(y): 
        if validregion[yii,xii]:
            nlm_ = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients
            nlm_[0], nlm_[3], nlm_[10] = 1, x_, y_
            m, t = np.array([0,0,1]), np.array([1,0,0])
            mm, mt = np.tensordot(m,m, axes=0), np.tensordot(m,t, axes=0)
            tau_ps_mm = 1*(np.identity(3)-3*mm) 
#            Ezz[yii,xii] = sf.Eeiej(nlm_, e1,e2,e3, Ecc,Eca,alpha,nprime)[-1,-1]
            Ezz[yii,xii] = sf.Evw(nlm_, mm,tau_ps_mm, Ecc,Eca,alpha,nprime)
        else:
            Ezz[yii,xii] = np.nan
                
print('done')

lvls = np.arange(0.2,1.6+1e-3,0.4)
CS = plt.contour(x,y,Ezz, levels=lvls, linewidths=0.7, linestyles=':', colors='0.2')
manual_locations = [(0.51, y_) for y_ in np.linspace(-0.5, 1.5, len(lvls))]
ax.clabel(CS, CS.levels, inline=True, fmt=r'$E_{zz}=%.1f$', fontsize=FS-3, manual=manual_locations)

### Model lines

hwmul=0.75
h_cdrx1 = plot_trajectory(ax, nlm_cdrx1, arrpos=15, c=c_cdrx, hwmul=hwmul)
h_cdrx2 = plot_trajectory(ax, nlm_cdrx2, arrpos=13, c=c_cdrx, hwmul=hwmul)

h_ddrx1 = plot_trajectory(ax, nlm_ddrx1, arrpos=40, c=c_ddrx, hwmul=hwmul)
h_ddrx2 = plot_trajectory(ax, nlm_ddrx2, arrpos=22, c=c_ddrx, hwmul=hwmul)
h_ddrx3 = plot_trajectory(ax, nlm_ddrx3, arrpos=28, c=c_ddrx, hwmul=hwmul)
h_ddrx4 = plot_trajectory(ax, nlm_ddrx4, arrpos=46, c=c_ddrx, hwmul=hwmul)

h_ss = plot_trajectory(ax, nlmr_ss, arrpos=None, c=c_smax, ls='--')
h_cc = plot_trajectory(ax, nlm_cc, arrpos=9,  c=c_smax)
h_ce = plot_trajectory(ax, nlm_ce, arrpos=17, c=c_girdle)

h_modellines = [h_ce, h_cc, h_ss, h_ddrx1, h_cdrx1]
legend_strings = ['Lattice rotation, unconf. extension', 'Lattice rotation, unconf. compression', 'Lattice rotation, simple shear', 'DDRX', 'CDRX']
legend_modellines = plt.legend(h_modellines, legend_strings, title=r'{\bf Modelled fabric state trajectories}', title_fontsize=FSLEG, loc=2, ncol=1, fontsize=FSLEG, frameon=False, **legkwargs)

### End-member cases

# Isotropic state
ax.plot(0,0,'o', ms=mse, c='k', label=None, zorder=20)
dytext = 0.04/normfac
plt.text(0, 0+dytext, r'{\bf Isotropic}', color='k', ha='center', va='bottom', fontsize=FSANNO)

# Unidirectional/delta-function (single max)
n20_delta = np.real(sp.sph_harm(0, 2, 0,0))/normfac
n40_delta = np.real(sp.sph_harm(0, 4, 0,0))/normfac
ax.plot(n20_delta,n40_delta, marker='o', ms=mse, ls='none', c=c_smax, label=None)
plt.text(n20_delta, n40_delta+dytext, '{\\bf Unidirectional}', color=c_smax, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# Planar (great circl) isotropy 
x, y = np.array([1,0,0]), np.array([0,1,0])
x2, y2 = np.einsum('i,j',x,x), np.einsum('i,j',y,y)
x4, y4 = np.einsum('i,j,k,l',x,x,x,x), np.einsum('i,j,k,l',y,y,y,y)
xy_sym = np.einsum('i,j',x,y) + np.einsum('i,j',y,x)
xy2 = x2 + y2
a2 = x2/2 + y2/2
a4 = x4/4 + y4/4 + np.einsum('ij,kl',xy2,xy2)/8 + np.einsum('ij,kl',xy_sym,xy_sym)/8 
nlm_girdle = np.real(sf.a4_to_nlm(a4))
#print(nlm_girdle, nlm_girdle[3],nlm_girdle[10])
x_, y_ = np.real(nlm_girdle[3])/normfac, np.real(nlm_girdle[10])/normfac
ax.plot(x_, y_, marker='o', ms=mse, ls='none', c=c_girdle, label=None)
plt.text(x_, y_+dytext, '{\\bf Planar}\n\\bf{isotropic}', color=c_girdle, ha='center', va='bottom', ma='center', fontsize=FSANNO)

# DDRX steady state
x_, y_ = np.real(nlm_ddrx2[-1,3]),np.real(nlm_ddrx2[-1,10])
ax.plot(x_, y_, marker='o', ms=mse, fillstyle='full', ls='none', c=c_ddrx, label=None)
plt.text(x_+0.01, y_-dytext, '{\\bf DDRX}\n\\bf{steady state}', color=c_ddrx, ha='center', va='top', ma='center', fontsize=FSANNO)

# Shading labels
plt.text(0.3/normfac, 0.18/normfac, '{\\bf Single maximum}', color=c_smax, ha='center', rotation=37, fontsize=FSANNO)
plt.text(-0.175/normfac, 0.125/normfac, '{\\bf Girdle}', color=c_girdle, ha='center', rotation=-40, fontsize=FSANNO)
plt.text(0.025/normfac, -0.325/normfac, '{\\bf Circle}', color=c_ddrx, ha='center', rotation=-25, fontsize=FSANNO)

plt.text(-0.25/normfac, -0.35/normfac, '{\\bf Unphysical}\n\\bf{eigenvalues}', color='0.3', ha='center', rotation=0, fontsize=FSANNO)

### Plot ODF insets?

if 1:

    inclination = 50 # view angle
    rot0 = -90
    rot = -20 + rot0 
    prj = ccrs.Orthographic(rot, 90-inclination)
    geo = ccrs.Geodetic()     

    W = 0.13 # ax width
    tickintvl=1

    ### Load measured data
    
    df = pd.read_csv('observed-states/Priestley-007.ctf.csv')
    sphcoords = qt.as_spherical_coords(qt.as_quat_array(df.to_numpy()[:,:-1]))
#    sphcoords = qt.as_spherical_coords(qt.as_quat_array(df.to_numpy()[::30,:]))
    qcolat, qlon = sphcoords[:,0], sphcoords[:,1] # COLAT [0;pi], LON [0;2pi]
    # Add reflected vector, too, for consistency
    qcolat = np.hstack((qcolat, np.pi-qcolat))
    qlon   = np.hstack((qlon, qlon-np.pi))
    qlat = np.pi/2 - qcolat 
    # Spectral coefs
    lm, nlm_len = sf.init(4)
    nlm_L4 = np.zeros((nlm_len), dtype=np.complex64) # The expansion coefficients
    caxes = np.array([ [np.cos(p)*np.sin(t), np.sin(p)*np.sin(t), np.cos(t)] for t, p in zip(qcolat,qlon) ])
    a2 = np.array([ np.einsum('i,j',c,c)         for c in caxes]).mean(axis=0)
    a4 = np.array([ np.einsum('i,j,k,l',c,c,c,c) for c in caxes]).mean(axis=0)
    nlm_L4[:16] = sf.a4_to_nlm(a4)
    nlm_L4 /= normfac
    # Rotated frame 
    (v1_colat, v1_lon, _) = get_v_angles(a2)
    nlmr_L4 = nlm_L4
    nlmr_L4 = sf.rotate_nlm(nlm_L4, 0, -v1_lon)
    nlmr_L4 = sf.rotate_nlm(nlmr_L4, +v1_colat, 0)
    # Rotated caxes
    Ry = np.matrix([[np.cos(-v1_colat),0,np.sin(-v1_colat)],[0,1,0],[-np.sin(-v1_colat),0,np.cos(-v1_colat)]]) # R_y
    Rz = np.matrix([[np.cos(-v1_lon),-np.sin(-v1_lon),0],[np.sin(-v1_lon),np.cos(-v1_lon),0],[0,0,1]]) # R_z
    caxes = np.array([np.einsum('ij,jk,k',Ry,Rz,cax) for cax in caxes])
    qlonr   = np.arctan2(caxes[:,1],caxes[:,0]) # arctan2(y,x)
    qcolatr = np.arccos(caxes[:,2])  
    qlatr   =  np.pi/2 - qcolatr
    
    ### Modeled data
    
    arrmag = 0.1
    arr = lambda ang: arrmag*np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
    ODF_plots = (\
        {'nlm':nlm_cc[int(Nt*5/10),:],    'title':'Model',    'cax':None, 'axloc':(0.455, 0.67), 'darr':arr(180), 'lvlmax':0.8}, \
        {'nlm':nlm_ce[int(Nt*6/10),:],    'title':'Model',    'cax':None, 'axloc':(0.19, 0.51), 'darr':arr(30),  'lvlmax':0.5}, \
        {'nlm':nlm_ddrx4[int(Nt*4/10),:], 'title':'Model',    'cax':None, 'axloc':(0.45, 0.16), 'darr':arr(0),   'lvlmax':0.45}, \
        {'nlm':nlmr_L4,                   'title':'Priestley','cax':(qlatr,qlonr,expr_Priestley['color']), 'darr':arr(-50), 'axloc':(0.60, 0.43), 'lvlmax':0.8}, \
    )

    for ODF in ODF_plots:

        axpos = [ODF['axloc'][0],ODF['axloc'][1], W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        nlm = ODF['nlm']*normfac
        lvls = np.linspace(0.0,ODF['lvlmax'],6)
    
        F, lon,lat = discretize_ODF(nlm, lm)
        F[F<0] = 0 # fix numerical/truncation errors
        h = axin.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend=('max' if lvls[0]==0.0 else 'both'), cmap='Greys', nchunk=5) # "nchunk" argument must be larger than 0 for constant-ODF (e.g. isotropy) is plotted correctly.

        # Arrow to ODF state
        n20_, n40_ = np.real(nlm[3])/normfac, np.real(nlm[10])/normfac
        ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                        xytext=(n20_+ODF['darr'][0]/normfac, n40_+sc**2*ODF['darr'][1]/normfac), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", facecolor='black'),zorder=20)            

        # Add grid lines
        kwargs_gridlines = {'ylocs':np.arange(-90,90+30,30), 'xlocs':np.arange(0,360+45,45), 'linewidth':0.5, 'color':'black', 'alpha':0.25, 'linestyle':'-'}
        gl = axin.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
        gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

        if ODF['cax'] is not None:
            qlatd, qlond = get_deg(ODF['cax'][0], ODF['cax'][1])
            axin.plot(qlond, qlatd, ls='none', marker='o', markersize=0.16, color=ODF['cax'][2], transform=geo) 
       
        axin.set_title(ODF['title'], fontsize=FS-3)

### Experimental data points

for ii,expr in enumerate(experiments):
    x = correlations[ii]['n20'] # n_2^0
    y = correlations[ii]['n40'] # n_4^0
    if expr['type']=='ss':  mrk = 's'
    if expr['type']=='ue':  mrk = '^'
    if expr['type']=='uc':  mrk = 'd'
    if expr['type']=='ucw': mrk = 'X'
    ax.plot(x,y, ms=ms, ls='none', color=expr['color'], fillstyle='none', marker=mrk, label=expr['plotname'], zorder=1+(len(experiments)-ii))

### Aux

plt.sca(ax)
plt.xlabel(r'$\hat{\psi}_2^0$')
plt.ylabel(r'$\hat{\psi}_4^0$')

leg = plt.legend(bbox_to_anchor=(1.41,1.02), title=r'{\bf Observations}', title_fontsize=FSLEG, fontsize=FSLEG, frameon=True, **legkwargs); 
leg.get_frame().set_linewidth(0.7);
ax.add_artist(legend_modellines)

### Limits

plt.xlim(xlims)
plt.ylim(ylims)
#xticks = np.arange(xlims[0],xlims[1]+1e-3,0.25)
#ax.set_xticks(xticks[2::2])
#ax.set_xticks(xticks[::1], minor=True)

### Second x axis for a_zz^(2) comparrison 

secax = ax.secondary_xaxis('top', functions=(n20_to_azz, azz_to_n20))
secax.set_xlabel(r'$a^{(2)}_{zz}$')
xticks = np.arange(0,1+1e-3,0.1)
secax.set_xticks(xticks[0::2])
secax.set_xticks(xticks[::1], minor=True)

### Save figure

#plt.tight_layout()
plt.savefig('%s.png'%(SELFNAME), dpi=175)
plt.close()

