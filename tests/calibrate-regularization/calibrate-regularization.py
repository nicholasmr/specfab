# N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

"""
Find the most aggressive hyper regularization exponent for a given L that ensures
single maxima, when developed under uniaxial compression up until strain = -0.99, have 
well-bounded power spectra (i.e. equal or less than the delta function spectrum)
"""

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
from numpy import linalg as LA
from scipy.optimize import minimize

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import discrete as sfdsc
from specfabpy import integrator as sfint
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex()

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.gridspec as gridspec

import warnings
warnings.filterwarnings("ignore")

#---------------------
# Run options
#---------------------

PLOT_SUMMARY = 1

SAVE_F90 = False

if len(sys.argv) == 2:
    L_list = [int(sys.argv[1]), ]
    
elif len(sys.argv) == 3:
    L_list = np.arange(int(sys.argv[1]), int(sys.argv[2])+1, 2)
    SAVE_F90 = True
    
else:
    L_list = [12,]    

#---------------------
# Setup
#---------------------

mod = dict(type='ps', axis=2, T=+1, r=0) # Mode of deformation
strain_target = -0.99 # calibration target

D, W = sf.ugrad_to_D_and_W(sf.pureshear_ugrad(mod['axis'],mod['r'],mod['T']))

Nt = 50 # calibration is target state this number of integration steps (should be relatively low, pushing it a bit!)

iota = +1 # deck of cards behaviour for lattice rotation

verboseint = 0 # verbose sfint integration

apply_bounds = False

#---------------------
# Initial guess, nu0
#---------------------

expo_opt = np.zeros((len(L_list)))
nu_opt   = np.zeros((len(L_list)))

for Lii, L in enumerate(L_list):

    if L == 4:  
        expo = 1.7
        nu0 = 1 

    if L == 6:  
        expo = 1.15
        nu0 = 1
       
    if L == 8:
        expo = 1.6
#        expo = 2.7 # aggressive, but seems to work? Maybe use this for FEM models
        nu0 = 1
        
    if L == 10:
        expo = 2
        nu0 = 5

    if L == 12:
        expo = 2
        nu0 = 5
        
    if L == 14:
        expo = 2
        nu0 = 10
        
    if L == 16:
        expo = 2.5
        nu0 = 10
        
    if L == 18:
        expo = 2.5
        nu0 = 10
        
    if L == 20: 
        expo = 3
        nu0 = 10
        
    #---------------------
    # Minimization problem, J(nu)
    #---------------------

    lm, nlm_len = sf.init(L)
    
    Sl_dirac, Lrange, nlm_dirac = sfdsc.Sl_delta(lm, sf)
    #print('*** S(l) for Delta function (l=0,2,...,L) = ', Sl_dirac)

    llrange = [int(ll*(ll+1)/2) for ll in Lrange] # (l,m=0) indices: 0, 3, 10, ....

    # require n20, n40 closely match the delta function solution.
    I = llrange[0:3] 
    
    def J(x, *args):
        nu = x[0]*nu0 
        nlm, *_ = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=iota, nu=nu, regexpo=expo, apply_bounds=apply_bounds, verbose=verboseint) 
        nlm_end = nlm[-1,:]
        err = np.real(nlm_dirac - nlm_end)[I] 
        abserr = np.abs(err)
        J = np.sum(abserr)
        print('x=%.8e :: nu=%.10e'%(x[0], nu))
        return J
    
    #---------------------
    # Estimate nu by solving inverse problem
    #---------------------
    
    print('*** Inverting for nu by minimizing the error in component with index = ', I)
    res = minimize(J, 1, method='L-BFGS-B', bounds=([0.1,10],), options=dict(eps=1e-4))

    nu = res.x[0]*nu0
    print(res)
    print('\n')
    print('*** Solution for L=%i is (nu0=%.3e):'%(nu0,L))
    print('  expo = %.3f'%(expo))
    print('  nu   = %.16f'%(nu))

    nu_opt[Lii]   = nu
    expo_opt[Lii] = expo
        
    #---------------------
    # Sample solution
    #---------------------
        
    nlm, F, time, *_ = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, iota=iota, nu=nu, regexpo=expo, apply_bounds=apply_bounds, verbose=verboseint) 
    dt = time[1]-time[0]
    tsteps = range(Nt+1)    
    strain_t = np.array([ sf.F_to_strain(sf.pureshear_F(mod['axis'],mod['r'],mod['T'], nn*dt))[2,2] for nn in tsteps])
    
    #---------------------
    # Plot results
    #---------------------
    
    if PLOT_SUMMARY:
        
        scale=0.825
        fig = plt.figure(figsize=(11*scale,7*scale))
        gs = gridspec.GridSpec(2, 3, height_ratios=[1,0.6])
        
        #---------------------
        
        ax = fig.add_subplot(gs[0,0])
        for nn in tsteps:
            Sl_model = np.array([sf.Sl(nlm[nn,:], l) for l in Lrange]) 
            Sl_model /= Sl_model[0] # normalize
            h = ax.semilogy(Lrange, Sl_model)
            
        ax.semilogy(Lrange, Sl_dirac, '--k', lw=1.5, label='$n=\delta$')  
        
        ax.set_xlabel('$l$')
        ax.set_ylabel('$S(l)/S(0)$')
        ax.set_ylim([1e-1,1.5])
        ax.set_xticks(np.arange(0,21,2))
        ax.set_xlim([0,np.amax([12,L])])
        ax.grid()
        ax.set_title('$L=%i$, $N_t=%i$'%(L, Nt))
        ax.legend()
        
        #---------------------
        
        vecdim = (3,Nt+1)
        eigvals = np.zeros(vecdim, dtype=np.float64)
        Eeiej = np.zeros((6,Nt+1), dtype=np.float64)
        for nn in tsteps:
            c = nlm[nn,:]
            e1,e2,e3, eigvals[:,nn] = sf.frame(c, 'e')
            (Eij_grain, alpha, n_grain) = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
            try: Eeiej[:,nn] = np.transpose(sf.Eij_tranisotropic(c, e1,e2,e3, Eij_grain, alpha, n_grain))
            except: pass
        
        ax = fig.add_subplot(gs[0,2])
        ax.plot(strain_t, eigvals[0,:], '-',  c='tab:red',   lw=2, label=r'$\lambda_1$')
        ax.plot(strain_t, eigvals[1,:], '-',  c='tab:blue',  lw=2, label=r'$\lambda_2$')
        ax.plot(strain_t, eigvals[2,:], '--', c='tab:green', lw=2, label=r'$\lambda_3$')
        ax.set_ylabel(r'$\lambda_i$')
        ax.set_xlabel(r'$\epsilon_{zz}$')
        ax.set_ylim([0,1])
        xlims_strain = [0,-1]
        ax.set_xlim(xlims_strain)
        ax.legend()
        ax.grid()
        ax.set_title('Eigenvalues')
        ev = eigvals[:,-1]
        print('*** Final eigenvalues are (sum %.3f): '%(sum(ev)), ev)
        
        #---------------------
        
        ax = fig.add_subplot(gs[0,1])
        ax.semilogy(strain_t, Eeiej[0], '-',  c='tab:red',   label=r'$E_{11}$')
        ax.semilogy(strain_t, Eeiej[2], '--', c='tab:blue',  label=r'$E_{33}$')
        ax.semilogy(strain_t, Eeiej[4], '-',  c='tab:green', label=r'$E_{13}$')
        ax.set_ylim([1e-2,1e1])
        ax.set_xlim(xlims_strain)
        ax.set_ylabel('$E_{ij}$')
        ax.set_xlabel('$\epsilon_{zz}$')
        ax.legend()
        ax.grid()
        ax.set_title('Eigenenhancements')
        
        #---------------------
        
        geo, prj = sfplt.getprojection(rotation=45, inclination=50)
        ax_ODF = [fig.add_subplot(gs[1, ii], projection=prj) for ii in [0,1,2]]
        
#        for ii, nn in enumerate([int(Nt/3), int(0.65*Nt), Nt]):
        if 1:
            ii,nn=2,Nt
            ax_ODF[ii].set_global()
            sfplt.plotODF(nlm[nn,:], lm, ax_ODF[ii], lvlset=(np.linspace(0,0.4,9),lambda x,p:'%.1f'%x), cblabel=r'$n(\epsilon_{zz}=%.2f)/N$'%(strain_t[nn]))
            sfplt.plotcoordaxes(ax_ODF[ii], geo, axislabels='vuxi')
            
        #---------------------
        
        gs.tight_layout(fig)
        fname = 'calibrate-regularization-L%i.png'%(L)
        print('*** Saving summary %s'%(fname))
        plt.savefig(fname, dpi=250)
    
    
if SAVE_F90:

    print('*** Dumping results to .f90 include file, used by dynamics.f90')
    
    with open('regcalib.f90', 'w') as f:
        for Lii, L in enumerate(L_list):
            f.write('if (Lcap == %i) then\n'%(L))
            f.write('    expo = %.3fd0\n'%(expo_opt[Lii]))
            f.write('    nu   = %.16fd0\n'%(nu_opt[Lii]))
            f.write('end if \n')
        
