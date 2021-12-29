# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
from numpy import linalg as LA
from scipy.optimize import minimize

sys.path.insert(0, '../../demo')
from header import *
from specfabpy import specfabpy as sf

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc

import warnings
warnings.filterwarnings("ignore")

#get_ipython().magic('clear')

#---------------------
# Setup
#---------------------

Nc = 50 # Number of "integration time steps" (we consider the analytical solution here).

TEST_GIRDLE = 0 # Validate calibration (for single maximum fabrics) against girdle fabrics.

if not TEST_GIRDLE: 
    epszz = -0.98 # Target: 2% of initial parcel height
    D = np.diag([0.5,0.5,-1]); Dzz = D[2,2];
    
else:
    epszz = 10 # Target: 10 times taller parcel compared to initial parcel height
    D = -np.diag([0.5,0.5,-1]); Dzz = D[2,2];

te = 1/Dzz # characteristic time, t_e
t_epszz = te*np.log(epszz + 1)
dt = t_epszz/Nc # time step size for thresshold strain epszz in "Nc" time steps 

tsteps = np.arange(Nc+1) # list of steps
epszz_t = np.exp(tsteps*dt/te)-1 # list of strains
print(Nc,dt)

W = 0*D # Assume strain-rate tensor is in eigen frame

#---------------------
# Initial guess, nu0
#---------------------

SOLVE_FOR_NU = 1 and (not TEST_GIRDLE)

L_list = [4,6,8]
L_list = [8,]
    
for L in L_list:

    nu0 = 2e0 # init guess
    
    if L == 4: expo = 1.65
    if L == 6: expo = 2.35
    if L == 8: expo = 3.0

    
    if TEST_GIRDLE:
        if L == 4: nu0 = 2.126845e+00
        if L == 6: nu0 = 2.892680e+00
        if L == 8: nu0 = 3.797282e+00
        
    #---------------------
    # Init specfab
    #---------------------
    
    lm, nlm_len  = sf.init(L)
    nlm      = np.zeros((nlm_len, Nc+1), dtype=np.complex128)
    nlm[0,:] = 1/np.sqrt(4*np.pi) # Init with isotropy
    nlm0 = nlm[:,0].copy()
    
    M = sf.dndt_LATROT(nlm0, D, W) # strain-rate assumed constant
    if 0 and TEST_GIRDLE: # test specfab implemention of solution.
        print("*** USING SPECFAB'S dndt_REG() (i.e. verifying implementation)")
        nu0 = 1
        R = sf.dndt_REG(nlm0, D) 
    else:
        L_ = sf.Lmat(nlm0)/(L*(L+1)) # normalized laplacian matrix
        R = -LA.norm(D) * np.power(np.abs(L_), expo)
    
    #print(np.diag(nu0*R))
    
    #---------------------
    # Solution
    #---------------------
    
    def f_nlm(nu_, Nc_):
        if 0: # Analytical solution
            w, V = LA.eig(M + nu_*R)
            t = dt*Nc_
            H = np.matmul(V,np.matmul(np.diag(np.exp(t * w)), LA.inv(V)))
            return np.matmul(H, nlm0), w
        else: # Numerical (Euler) solution
            nlm_prev = nlm0.copy()
            for ii in np.arange(Nc_):
                nlm_next = nlm_prev + dt*np.matmul(M + nu_*R,nlm_prev)
                nlm_prev = nlm_next.copy()
            return nlm_next, 0
        
    #---------------------
    # Power spectrum, S(l)
    #---------------------
    
    Lrange = np.arange(0,L+1,2) # 0, 2, 4, 6, ...
    llrange = [int(ll*(ll+1)/2) for ll in Lrange] # indices 0, 3, 10, ....
    powerspectrum = lambda l,cl0: 1/(2*l+1) * np.abs(cl0)**2
    Sl0 = powerspectrum(0, nlm0[0])   
    def Sl(nlm): return 1/Sl0 * np.array([ powerspectrum(L,nlm[llrange[ii]]) for ii, L in enumerate(Lrange) ]) 
    
    nlm_dirac = nlm0.copy()
    for ii, lli in enumerate(llrange):
        nlm_dirac[lli] = sp.sph_harm(0, Lrange[ii], 0,0)
    Sl_dirac = Sl(nlm_dirac)
    
    #---------------------
    # Minimization problem, J(nu)
    #---------------------
    
    I = llrange[0:2] # require nlm entries 1 (=n20), 2 (=n40) closely match the delta function solution.
    #print(Lrange[I])
    
    def J(nu_test, *args):
        (Nc_, ) = args
        nlm, _ = f_nlm(nu_test[0], Nc_)
        err = np.abs(nlm-nlm_dirac)[I]
        return np.sum(err) # = J
    
    #---------------------
    # Solve problem
    #---------------------
    
    if not SOLVE_FOR_NU:
        nu = nu0    
        print('*** SPECIFIED nu0, nu = ', nu0, nu)
        
    else:
        res = minimize(J, nu0, method='L-BFGS-B', args=(Nc), options={'disp': False})
        nu = res.x[0]
        print(res)
        print('*** SOLVED (nu0=%.3e):'%(nu0))
        print('\texpo = %.2f'%(expo))
        print('\tnu   = %e'%(nu))
        
    #---------------------
    # Sample solution
    #---------------------
        
    for nn in np.arange(Nc): nlm[:,1+nn], _ = f_nlm(nu, nn+1)
    
    #---------------------
    # Plot results
    #---------------------
    
    scale=0.9
    fig = plt.figure(figsize=(11*scale,7*scale))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1,0.7])
    
    #---------------------
    
    ax = fig.add_subplot(gs[0,0])
    for nn in tsteps:
        h = ax.semilogy(Lrange, Sl(nlm[:,nn]))
        
    ax.semilogy(Lrange, Sl_dirac, '--k', lw=1.5, label='$\delta(\\vu{r}-\\vu{z})$')  
    
    ax.set_xlabel('$l$')
    ax.set_ylabel('$S(l)/S(0)$')
    ax.set_ylim([1e-1,1.5])
    ax.set_xticks(np.arange(0,13,2))
    ax.set_xlim([0,12])
    ax.grid()
    ax.set_title('$L=%i$, $N_c=%i$'%(L, Nc))
    
    #---------------------
    
    vecdim = (3,Nc+1)
    eigvals = np.zeros(vecdim, dtype=np.float64)
    Eeiej = np.zeros((3,3,Nc+1), dtype=np.float64)
    for nn in tsteps:
        c = nlm[:,nn]
        e1,e2,e3, eigvals[:,nn] = sf.frame(c, 'e')
        #print(nn, eigvals[:,nn])
        Eca_lin,  Ecc_lin,  alpha_lin, nprime  = sf.Eca_opt_lin,  sf.Ecc_opt_lin,  sf.alpha_opt_lin, 1  # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
        Eeiej[:,:,nn] = np.transpose(sf.Eeiej(c, e1,e2,e3, Ecc_lin, Eca_lin, alpha_lin, nprime))
    
    ax = fig.add_subplot(gs[0,2])
    ax.plot(epszz_t, eigvals[0,:])
    ax.plot(epszz_t, eigvals[1,:])
    ax.plot(epszz_t, eigvals[2,:])
    ax.set_ylabel('$a_i$')
    ax.set_xlabel('$\epsilon_{zz}$')
    ax.set_ylim([0,1])
    ax.set_xlim(epszz_t[[0,-1]])
    ax.grid()
    
    #---------------------
    
    ax = fig.add_subplot(gs[0,1])
    ax.semilogy(epszz_t, Eeiej[-1,-1])
    ax.semilogy(epszz_t, Eeiej[0,-1])
    ax.set_ylim([1e-2,1e1])
    ax.set_xlim(epszz_t[[0,-1]])
    ax.set_ylabel('$E_{ij}$')
    ax.set_xlabel('$\epsilon_{zz}$')
    ax.grid()
    
    #---------------------
    
    inclination = 50 # view angle
    rot0 = -90
    rot = -40 + rot0 
    prj = ccrs.Orthographic(rot, 90-inclination)
    geo = ccrs.Geodetic()     
    ax_ODF = [fig.add_subplot(gs[1, ii], projection=prj) for ii in [0,1,2]]
    lm, nlm_len = sf.init(L)
    for ii, nn in enumerate([int(Nc/3), int(0.65*Nc), Nc]):
        epszz_nn = np.exp(dt*nn/te) - 1
        plot_ODF(nlm[:nlm_len, nn], lm, ax=ax_ODF[ii], cmap='Greys', cblabel=r'$\psi(\epsilon_{zz}=%.2f)$'%(epszz_nn), latres=40)
        ax_ODF[ii].set_global()
    
    #---------------------
    
    gs.tight_layout(fig)
    fname = '%s_L%i.png'%('girdle' if TEST_GIRDLE else 'smax', L)
    print('*** Saving summary %s'%(fname))
    plt.savefig(fname, dpi=150)