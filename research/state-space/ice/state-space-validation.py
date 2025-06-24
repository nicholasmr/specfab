# N. M. Rathmann <rathmann@nbi.ku.dk> and D. Lilien, 2022-2023

"""
CPO state-space map comparing model trajectories to data (Lilien et al., 2023)
"""

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import pandas as pd
import pickle, glob

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors

sys.path.insert(0, '..')
from localheader import *
from experiments import * # experiment definitions (data structures for experiment files, etc.)

from specfabpy import specfab as sf
from specfabpy import constants as sfconst
from specfabpy import integrator as sfint
from specfabpy import statespace as sfsp
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=12)
FSLEG = FS-1.5
FSANNO = FS-2.5

SELFNAME = sys.argv[0][:-3] # used as prefix for pickled files
os.system('mkdir -p specfab-state-trajectories')
def pfile(fname): return "specfab-state-trajectories/%s--%s.p"%(SELFNAME, fname) # full path name for pickled files

#--------------------
# Run options
#--------------------

DEBUG           = 0  # Low-resolution plotting, etc.
INTEGRATE_MODEL = 0  # Generate model lines? Else load saved Pickle files.
            
PLOT_ODF_INSETS = True
            
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
    
lm, nlm_len = sf.init(L) 
    
#--------------------
# Experimental data to plot
#--------------------

# All data (Gerber et al., 2023)
experiments = (expr_GRIP, expr_LAWDOME, expr_EGRIP_MAIN, expr_SPICE, expr_EGRIP_S5, expr_Priestley, expr_Qi, expr_Fan_10,expr_Fan_4, expr_Hunter) 

# Lilien et al. (2022)
#experiments = (expr_GRIP, expr_LAWDOME, expr_EGRIP_MAIN, expr_SPICE, expr_Priestley, expr_Qi, expr_Fan_10,expr_Fan_4, expr_Hunter) 

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
    corr['n20'] /= norm # Normalization
    corr['n40'] /= norm
    corr['n60'] /= norm
    correlations.append(corr) 
    
#--------------------
# Modeled correlations
#--------------------

def integrate_model(nlm0, Mtype, modtype, Nt=Nt, rotate=False, name=None):

    iota = zeta = Gamma0 = Lambda = nu = None

    if Mtype == 'LROT': iota, zeta, nu = 1, 0, 1
    if Mtype == 'DDRX': Gamma0 = 18e-0
    if Mtype == 'CDRX': Lambda = 2e-1
        
    if modtype == 'uc': mod, target = dict(type='ps', axis=2, T=+1, r=0), -0.99
    if modtype == 'ue': mod, target = dict(type='ps', axis=2, T=-1, r=0), 6
    if modtype == 'ss': mod, target = dict(type='ss', plane=1, T=+1), np.deg2rad(83)
    
    nlm, F, *_ = sfint.lagrangianparcel(sf, mod, target, Nt=Nt, nlm0=nlm0, iota=iota, zeta=zeta, Gamma0=Gamma0, Lambda=Lambda, nu=nu) 

    if rotate:
        nlmr = nlm.copy()
        for tt in np.arange(1,Nt+1):
            (m1_colat, m1_lon, _) = get_m_angles(sf.a2(nlm[tt,:]))
            nlmr[tt,:] = sf.rotate_nlm(sf.rotate_nlm( nlm[tt,:], 0, -m1_lon), -m1_colat, 0)
    else:
        nlmr = nlm # no rotation requested, just copy the unrotated solution

    if name is not None: 
        pickle.dump([nlm, nlmr, lm, nlm_len], open(pfile(name), "wb"))

    return nlm, nlmr    
    
    
if INTEGRATE_MODEL:

    print('*** Generating model trajectories from scratch. You can re-use the saves trajectories to avoid re-calculating them by setting INTEGRATE_MODEL=1')

    ### Solve for state-vector time evolution
    
    nlm_iso = np.zeros((nlm_len), dtype=np.complex64)
    nlm_iso[0] = norm
    
    nlm_uc, _  = integrate_model(nlm_iso, 'LROT', 'uc', name='uc') 
    nlm_ue, _  = integrate_model(nlm_iso, 'LROT', 'ue', name='ue') 
    _, nlmr_ss = integrate_model(nlm_iso, 'LROT', 'ss', name='ss', rotate=True) 

    f_nlm_init = lambda nlm, n20_ref: nlm[np.argmin(np.abs(nlm[:,sf.I20]-n20_ref)), :] 
    nlm_ddrx1, _ = integrate_model(f_nlm_init(nlm_ue, -0.20), 'DDRX', 'uc', name='ddrx1') # DDRX trajectory 1
    nlm_ddrx2, _ = integrate_model(f_nlm_init(nlm_uc, +0.00), 'DDRX', 'uc', name='ddrx2') # DDRX trajectory 2
    nlm_ddrx3, _ = integrate_model(f_nlm_init(nlm_uc, +0.20), 'DDRX', 'uc', name='ddrx3') # DDRX trajectory 3
    nlm_ddrx4, _ = integrate_model(f_nlm_init(nlm_uc, +0.40), 'DDRX', 'uc', name='ddrx4') # DDRX trajectory 4
    
    nlm_cdrx1, _ = integrate_model(f_nlm_init(nlm_ue, -0.27), 'CDRX', 'uc', name='cdrx1') # DDRX trajectory 1
    nlm_cdrx2, _ = integrate_model(f_nlm_init(nlm_uc, +0.32), 'CDRX', 'uc', name='cdrx2') # DDRX trajectory 1
    
    
### Load solutions

def load_solution(expname):
    nlm, nlmr,  lm, nlm_len = pickle.load(open(pfile(expname), "rb"))    
    nlm  = np.array([ nlm[tt,:]/nlm[tt,0]   for tt in np.arange(Nt+1) ]) # normalize
    nlmr = np.array([ nlmr[tt,:]/nlmr[tt,0] for tt in np.arange(Nt+1) ]) # normalize
    return (nlm, nlmr, lm, nlm_len)

nlm_uc, _,  lm, nlm_len = load_solution('uc')
nlm_ue, _,  lm, nlm_len = load_solution('ue')
_, nlmr_ss, lm, nlm_len = load_solution('ss')

nlm_ddrx1, _, lm, nlm_len = load_solution('ddrx1')
nlm_ddrx2, _, lm, nlm_len = load_solution('ddrx2')
nlm_ddrx3, _, lm, nlm_len = load_solution('ddrx3')
nlm_ddrx4, _, lm, nlm_len = load_solution('ddrx4')

nlm_cdrx1, _, lm, nlm_len = load_solution('cdrx1')
nlm_cdrx2, _, lm, nlm_len = load_solution('cdrx2')
                
#--------------------
# Construct plot
#--------------------

### Setup figure

scale = 4.7
fig = plt.figure(figsize=(1.6*scale,1.1*scale))
plt.subplots_adjust(left=0.08, right=0.725, top=0.90, bottom=0.125)
ax = plt.gca()

xlims, ylims = [-1.40,2.65], [-1.8,3.75]
sc = np.diff(ylims)/np.diff(xlims)

ms = 6.0
c_ddrx    = '#b2182b'
cl_circle = '#fddbc7'
c_cdrx   = 'k'

legkwargs = {'handlelength':1.4, 'framealpha':1.0, 'fancybox':False, 'handletextpad':0.4, 'borderpad':0.37, 'edgecolor':'k'}

### Determine valid subspace

grain_params = sfconst.ice['viscoplastic']['linear'] # Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
(Ezz, Exz, x, y, isvalid) = sfsp.Eij_statemap_ice(grain_params, xlims, ylims, resx=RESX, resy=RESY)

C = statespace_shading(nlm_uc, nlm_ue, x, y, isvalid, shade_circle=True)

# Plot valid subspace and fabric-type shadings 
im = ax.imshow(C, aspect='auto', extent=[np.amin(xlims), np.amax(xlims), np.amin(ylims), np.amax(ylims)], origin='lower', zorder=1)

# Plot E_zz
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

h_ss = plot_trajectory(ax, nlmr_ss, arrpos=None, c=sfsp.c_unidir, ls='--')
h_uc = plot_trajectory(ax, nlm_uc, arrpos=9,  c=sfsp.c_unidir)
h_ue = plot_trajectory(ax, nlm_ue, arrpos=17, c=sfsp.c_planar)

h_modellines = [h_ue, h_uc, h_ss, h_ddrx1, h_cdrx1]
legend_strings = ['Lattice rotation, uniaxial extension', 'Lattice rotation, uniaxial compression', 'Lattice rotation, simple shear', 'DDRX', 'CDRX']
legend_modellines = plt.legend(h_modellines, legend_strings, title=r'{\bf Modelled fabric state trajectories}', title_fontsize=FSLEG, loc=2, ncol=1, fontsize=FSLEG, frameon=False, **legkwargs)

### Labels

dy0 = 0.065
mse = 7.5
sfsp.plot_nlm_cases(ax, FSANNO, ms=mse, show_circle=False, dy0=dy0, isolbl_above=True)

nl0_unidir, nl0_planar, nl0_circle = sfdsc.nlm_ideal_cases(norm=norm)
ax.plot(nl0_circle[0], nl0_circle[1] , c=c_ddrx, marker='o', ms=mse, ls='none', label=None, zorder=20)
plt.text(nl0_circle[0], nl0_circle[1]-1.4*dy0/norm, '{\\bf DDRX}\n{\\bf steady state}', color=c_ddrx, ha='center', va='center', ma='center', fontsize=FSANNO)

kwargs_lbl = dict(ha='center', fontsize=FSANNO)
plt.text(+0.300/norm, +0.180/norm, '{\\bf Single maximum}', color=sfsp.c_unidir, rotation=37, **kwargs_lbl)
plt.text(-0.175/norm, +0.125/norm, '{\\bf Girdle}', color=sfsp.c_planar, rotation=-40, **kwargs_lbl)
plt.text(+0.025/norm, -0.325/norm, '{\\bf Circle}', color=c_ddrx, rotation=-25, **kwargs_lbl)
plt.text(-0.250/norm, -0.350/norm, '{\\bf Unphysical}\n{\\bf eigenvalues}', color=sfplt.c_dgray, rotation=0, **kwargs_lbl)

### Plot ODF insets?

if PLOT_ODF_INSETS:

    ### Define observational states
    
    (q, qr, m, mnew, caxes, nlm, nlmr, lm) = load_sample('007.ctf.csv', expr_Priestley)
    (qlatr, qcolatr, qlonr) = qr
    nlmr_L4 = nlmr[:sf.L4len]/norm
    cax_Priestley = (qlatr,qlonr,expr_Priestley['color'])
    
    ### Define model states
    
    arr = lambda ang: 0.1 * np.array([np.cos(np.deg2rad(ang)),np.sin(np.deg2rad(ang))])
    
    ODF_plots = (\
        {'nlm':nlm_uc[int(Nt*4.7/10),:],    'title':'Model',    'cax':None,          'axloc':(0.46, 0.675), 'darr':arr(180), 'lvlmax':0.8}, \
        {'nlm':nlm_ue[int(Nt*6/10),:],      'title':'Model',    'cax':None,          'axloc':(0.185, 0.51), 'darr':arr(30), 'lvlmax':0.5}, \
        {'nlm':nlm_ddrx4[int(Nt*3.5/10),:], 'title':'Model',    'cax':None,          'axloc':(0.45, 0.16), 'darr':arr(0), 'lvlmax':0.45}, \
        {'nlm':nlmr_L4,                     'title':'Priestley','cax':cax_Priestley, 'axloc':(0.60, 0.43), 'darr':arr(-50),'lvlmax':0.8}, \
    )

    ### Plot
    
    for ODF in ODF_plots:

        # Setup axis
        geo, prj = sfplt.getprojection(rotation=55+180, inclination=50)   
        W = 0.12 # axis width
        axin = plt.axes([ODF['axloc'][0],ODF['axloc'][1], W,W], projection=prj) #, transform=ax.transData)
        axin.set_global()
        axin.set_title(ODF['title'], fontsize=FS-3)
                
        # Plot ODF
        nlm = ODF['nlm']*norm
        sfplt.plotODF(nlm, lm, axin, lvlset=(np.linspace(0.0,ODF['lvlmax'],6), lambda x,p:'%.1f'%x), showcb=False)

        # Arrow to ODF state
        n20_ = np.real(nlm[sf.I20])/norm
        n40_ = np.real(nlm[sf.I40])/norm
        ax.annotate("", xy=(n20_, n40_), xycoords='data', \
                        xytext=(n20_+ODF['darr'][0]/norm, n40_+sc**2*ODF['darr'][1]/norm), textcoords='data', \
                        arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", facecolor='black'),zorder=20)            

        # Plot c-axis data
        if ODF['cax'] is not None:
            qlatd, qlond = np.rad2deg(ODF['cax'][0]), np.rad2deg(ODF['cax'][1])
            axin.plot(qlond, qlatd, ls='none', marker='o', markersize=0.16, color=ODF['cax'][2], transform=geo) 
       

### Experimental data points

for ii,expr in enumerate(experiments):
    x = correlations[ii]['n20'] # n_2^0
    y = correlations[ii]['n40'] # n_4^0
    if expr['type']=='ss':  mrk = 's'
    if expr['type']=='ue':  mrk = '^'
    if expr['type']=='uc':  mrk = 'd'
    if expr['type']=='ucw': mrk = 'X'
    ax.plot(x,y, ms=ms, ls='none', color=expr['color'], fillstyle='none', marker=mrk, label=expr['plotname'], zorder=1+(len(experiments)-ii))

### Axes and legends

plt.sca(ax)

plt.xlabel(r'$\hat{n}_2^0$')
plt.ylabel(r'$\hat{n}_4^0$')

leg = plt.legend(bbox_to_anchor=(1.43,1.02), title=r'{\bf Observations}', title_fontsize=FSLEG, fontsize=FSLEG, frameon=True, **legkwargs); 
leg.get_frame().set_linewidth(0.7);
ax.add_artist(legend_modellines)

xt = np.arange(-5,5,0.25)
ax.set_xticks(xt[::4])
ax.set_xticks(xt, minor=True)
ax.set_yticks(xt[::4])
ax.set_yticks(xt, minor=True)

plt.xlim(xlims)
plt.ylim(ylims)

### Second x axis for a_zz^(2) comparrison 

secax = ax.secondary_xaxis('top', functions=(n20_to_azz, azz_to_n20))
secax.set_xlabel(r'$a^{(2)}_{zz}$')
xticks = np.arange(0,1+1e-3,0.1)
secax.set_xticks(xticks[0::2])
secax.set_xticks(xticks[::1], minor=True)

### Save figure

plt.savefig('%s.png'%(SELFNAME), dpi=175)

