# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2022

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

#---------------------
# Setup
#---------------------

### Settings
SOLVE_FOR_NU = 1  # Run minimization to determine regularization magnitude per L
APPLY_BOUNDS = 0  # Test effect of applying bounds to solution (after solving for nu without bounds enabled)
TEST_GIRDLE  = 0  # Validate calibration (for single maximum fabrics) against girdle fabrics.

L_list = [4,6,8,20]
L_list = [8,10,12,]
#L_list = [10,]

### Mode of deformation
if not TEST_GIRDLE: 
    epszz = -0.99 # Target: 1% of initial parcel height
    D = np.diag([0.5,0.5,-1]);
    
else:
    epszz = 10 # Target: 10 times taller parcel compared to initial parcel height
    D = -np.diag([0.5,0.5,-1]); 

Dzz = D[2,2]; # size of compressive strain-rate
W = 0*D # Assume strain-rate tensor is in eigen frame

### Rotate compressive axis
latrot = np.deg2rad(0)
print('Strain-rate tensor sym. axis has orientation theta=%.4f deg.'%(np.rad2deg(latrot)))
c, s = np.cos(latrot), np.sin(latrot)
R = np.matrix([[c,0,s],[0,1,0],[-s,0,c]])
D = np.matmul(np.matmul(R,D),R.T)

### Numerics
Nc = 50 # Number of "integration time steps" (we consider the analytical solution here).
te = 1/Dzz # characteristic time, t_e
t_epszz = te*np.log(epszz + 1)
dt = t_epszz/Nc # time step size for thresshold strain epszz in "Nc" time steps 
tsteps = np.arange(Nc+1) # list of steps
epszz_t = np.exp(tsteps*dt/te)-1 # list of strains
print('Nc=%i, dt=%.4e s'%(Nc,dt))

#---------------------
# Initial guess, nu0
#---------------------

SOLVE_FOR_NU = SOLVE_FOR_NU and (not TEST_GIRDLE) # don't solve for nu for other modes of deformation than unconfined uniaxial compression.
    
for L in L_list:

    if L == 4:  
#        expo = 1.5
#        nu0 = 1 # init guess for iteration

        SOLVE_FOR_NU = 0
        expo = 1.55
        nu0 = 2.9 # init guess for iteration

    if L == 6:  
        #expo = 2.0
        #nu0 = 1 # init guess for iteration

        SOLVE_FOR_NU = 0
        expo = 1.25
        nu0 = 3.6
       
    if L == 8:
#        expo = 3 # 2.8
#        nu0 = 4e+0 # init guess for iteration
        expo = 1.5 # 2.8
        nu0 = 6.5e+0 # init guess for iteration
        
    if L == 10:
#        expo = 1.6
#        nu0 = 8e+0 # init guess for iteration
        expo = 1.7
        nu0 = 9e+0 # init guess for iteration
        
    if L == 12:
        expo = 1.7
        nu0 = 9e+0 # init guess for iteration
#        expo = 1.6
#        nu0 = 9e+0 # init guess for iteration
        
    if L == 20: 
        expo = 3
        nu0 = 3e+0 # init guess for iteration
    
    if TEST_GIRDLE:
        if L == 4: nu0 = 2.126845e+00
        if L == 6: nu0 = 2.892680e+00
        if L == 8: nu0 = 3.797282e+00
        
    #---------------------
    # Init specfab
    #---------------------
    
    lm, nlm_len = sf.init(L)
    nlm = np.zeros((nlm_len, Nc+1), dtype=np.complex128)
    nlm[0,:] = 1/np.sqrt(4*np.pi) # Init with isotropy
    nlm0 = nlm[:,0].copy()
    
    M = sf.M_LROT(nlm0, D, W) # strain-rate assumed constant for calibration experiments
    if 0 and TEST_GIRDLE: # test specfab implemention of solution.
        print("*** USING SPECFAB'S M_REG() (verifying fortran implementation)")
        nu0 = 1
        R = sf.M_REG(nlm0, D) 
    else:
        L_ = sf.Lmat(nlm0)/(L*(L+1)) # normalized laplacian matrix
        R = -LA.norm(D) * np.power(np.abs(L_), expo)
    
    #print(np.diag(nu0*R))
    
    #---------------------
    # Solution
    #---------------------

    def f_nlm(nu_, Nc_, apply_bounds=0):
    
        if 0: # Analytical solution (matrix exponentiation)
        
            w, V = LA.eig(M + nu_*R)
            t = dt*Nc_
            H = np.matmul(V,np.matmul(np.diag(np.exp(t * w)), LA.inv(V)))
            
            return np.matmul(H, nlm0), w
            
        else: # Numerical (Euler) solution
        
            nlm_prev = nlm0.copy()
            for ii in np.arange(Nc_):
                nlm_next = nlm_prev + dt*np.matmul(M + nu_*R,nlm_prev)
                if apply_bounds: nlm_next = sf.apply_bounds(nlm_next)
                nlm_prev = nlm_next.copy()
                
            return nlm_next, 0
        
    #---------------------
    # Delta function power spectrum, S_dirac(l)
    #---------------------
    
    nlm_dirac = nlm0.copy()
    for ii, (l,m) in enumerate(lm.T): nlm_dirac[ii] = sp.sph_harm(m,l, 0,latrot)
    Lrange = np.arange(0,L+1,2) # 0, 2, 4, 6, ...
    Sl_dirac = np.array([sf.Sl(nlm_dirac, l) for l in Lrange])
    Sl_dirac /= Sl_dirac[0] # normalize
    print('S_dirac(l) (for l=0,2,...,L) = ', Sl_dirac)
    
    #---------------------
    # Minimization problem, J(nu)
    #---------------------

    llrange = [int(ll*(ll+1)/2) for ll in Lrange] # indices 0, 3, 10, ....    
    I = llrange[1:2] # require nlm entries 1 (=n20), 2 (=n40) closely match the delta function solution.
    tol = 0.99
    
    def J(nu_test, *args):
        (Nc_, ) = args
        nlm, _ = f_nlm(nu_test[0], Nc_)
        abserr = np.abs( np.real(nlm - tol*nlm_dirac)[I] )
        return np.sum(abserr)
    
    #---------------------
    # Estimate nu by solving inverse problem
    #---------------------
    
    if not SOLVE_FOR_NU:
        nu = nu0    
        print('*** TESTING GUESS nu0, nu = ', nu0, nu)
        
    else:
        print('*** INVERTING for nu by minimizing the error in component with index = ', I, 'with tolerence = %.4f'%(tol))
        res = minimize(J, nu0, method='CG', args=(Nc), options={'disp': False})
        nu = res.x[0]
        print(res)
        print('*** SOLVED (nu0=%.3e):'%(nu0))
        print('        if (Lcap == %i) then'%(L))
        print('            expo = %.3f'%(expo))
        print('            nu   = %e'%(nu))
        
    #---------------------
    # Sample solution
    #---------------------
        
    for nn in np.arange(Nc): nlm[:,1+nn], _ = f_nlm(nu, nn+1, apply_bounds=APPLY_BOUNDS)
    
    #---------------------
    # Plot results
    #---------------------
    
    scale=0.825
    fig = plt.figure(figsize=(11*scale,7*scale))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1,0.6])
    
    #---------------------
    
    ax = fig.add_subplot(gs[0,0])
    for nn in tsteps:
        Sl_model = np.array([sf.Sl(nlm[:,nn], l) for l in Lrange]) 
        Sl_model /= Sl_model[0] # normalize
        h = ax.semilogy(Lrange, Sl_model)
        
    ax.semilogy(Lrange, Sl_dirac, '--k', lw=1.5, label='$\psi=\delta$')  
    
    ax.set_xlabel('$l$')
    ax.set_ylabel('$S(l)/S(0)$')
    ax.set_ylim([1e-1,1.5])
    ax.set_xticks(np.arange(0,21,2))
    ax.set_xlim([0,np.amax([10,L])])
    ax.grid()
    ax.set_title('Pwr spec. :: $L=%i$, $N_c=%i$'%(L, Nc))
    ax.legend()
    
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
    ax.plot(epszz_t, eigvals[0,:], '-',  c='tab:red',   lw=2, label=r'$\lambda_1$')
    ax.plot(epszz_t, eigvals[1,:], '-',  c='tab:blue',  lw=2, label=r'$\lambda_2$')
    ax.plot(epszz_t, eigvals[2,:], '--', c='tab:green', lw=2, label=r'$\lambda_3$')
    ax.set_ylabel(r'$\bf{a}^{(2)}$ eigen values')
    ax.set_xlabel(r'$\epsilon_{zz}$')
    ax.set_ylim([0,1])
    ax.set_xlim(epszz_t[[0,-1]])
    ax.legend()
    ax.grid()
    ax.set_title('Eigen values')
    print('eigvals for last step are: ', eigvals[:,-1])
    
    #---------------------
    
    ax = fig.add_subplot(gs[0,1])
    ax.semilogy(epszz_t, Eeiej[0,0], '-',  c='tab:red',   label=r'$E_{11}$')
    ax.semilogy(epszz_t, Eeiej[2,2], '--', c='tab:blue',  label=r'$E_{33}$')
    ax.semilogy(epszz_t, Eeiej[0,2], '-',  c='tab:green', label=r'$E_{13}$')
    ax.set_ylim([1e-2,1e1])
    ax.set_xlim(epszz_t[[0,-1]])
    ax.set_ylabel('$E_{ij}$')
    ax.set_xlabel('$\epsilon_{zz}$')
    ax.legend()
    ax.grid()
    ax.set_title('Eigen enhancements')
    
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
        if latrot==0: ax_ODF[ii].plot([0],[90], marker=r'$e_1$', ms=9, c='tab:orange', transform=geo) 
    #---------------------
    
    gs.tight_layout(fig)
    fname = 'calibrate-regularization-%s-L%i.png'%('girdle' if TEST_GIRDLE else 'singlemax', L)
    print('*** Saving summary %s'%(fname))
    plt.savefig(fname, dpi=250)
    
