# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020-2021

import sys, os, copy, code # code.interact(local=locals())
import numpy as np
import scipy.special as sp

sys.path.insert(0, '..')
from header import *
from specfabpy import specfabpy as sf # To use specfabpy compile the specfab Python module by running "make specfabpy"

import warnings
warnings.filterwarnings("ignore")
from progress.bar import Bar

year2sec = 31556926;

#----------------------
# Experiment
#----------------------

T_EXP_CC = 1 # confined vertical compression 
T_EXP_SS = 2 # vertical shear

T_EXP_STR = {T_EXP_CC:'cc', T_EXP_SS:'ss'}

### Select experiment

T_EXP = T_EXP_CC
T_EXP = T_EXP_SS

DEBUG = 0

#----------------------
# Parameters
#----------------------

nglen = 3
Aglen = 3.5e-26 # A(T=-25 deg.)
Gamma0 = 1*6e-6
if T_EXP==T_EXP_SS: DDRX_strainthres = 5.5
if T_EXP==T_EXP_CC: DDRX_strainthres = -0.7 


#----------------------
# Numerics
#----------------------

if T_EXP==T_EXP_SS:
    tau0mag = 1.087e6
    tau0 = tau0mag * np.matrix([[0,0,1],[0,0,0],[1,0,0]]) # Pa
    ugrad0 = 2*sf.eps_of_tau__isotropic(tau0, Aglen,nglen)[0,2]
    Nt = 2 * 200
                
if T_EXP==T_EXP_CC:
    tau0mag = 10e5
    r = 1
    tau0 = tau0mag * np.diag([-(1+r)/2., -(1-r)/2., 1]) # Pa
    ugrad0 = 2*sf.eps_of_tau__isotropic(tau0, Aglen,nglen)[2,2]
    Nt = 3 * 200
        
t_c = 1/ugrad0
T = 1.3 * year2sec
dt = T/Nt

if DEBUG:
    L = 12 # Spectral truncation
    nu0 = 3.0e-3 # Regularization strength (calibrated for L)
#    nu0 = 0.0e-14 # Regularization strength (calibrated for L)
else:
    L = 26 # Spectral truncation
    nu0 = 3e-3 # Regularization strength (calibrated for L)
        
        
#----------------------
# Grain parameters for enhancement-factor model
#----------------------

# Linear (n'=1) mixed Taylor--Sachs enhancements            
# Optimal n'=1 (lin) grain parameters (Rathmann and Lilien, 2021)
nprime    = 1 
Eca_lin   = sf.eca_opt_lin
Ecc_lin   = sf.ecc_opt_lin
alpha_lin = sf.alpha_opt_lin  

#----------------------
# Velocity gradient
#----------------------
    
if T_EXP==T_EXP_SS:
        
    SS = SimpleShear(t_c, plane='zx') 
    omg, eps = SS.W(), SS.D()
    print("*** eps_xz (for parcel deformation): ",eps[0,2])
    time = dt*np.arange(0,Nt)
    strain = SS.kappa * time
    ODF_strains = [0,2,5,7]
    
if T_EXP==T_EXP_CC:
        
    PS = PureShear(t_c, r, ax='z') 
    omg, eps = PS.W(), PS.D()
    print("*** eps_zz (for parcel deformation): ",eps[2,2])
    time = dt*np.arange(0,Nt)
    strain = np.array([PS.strain(tt)[2,2] for tt in time])
    ODF_strains = [0, -0.3, -0.6, -0.9] # maybe use DDRX_strainthres as one of the steps

ODF_tsteps = [np.argmin(np.abs(strain-thres)) for thres in ODF_strains]

print("*** strain max: ", strain[-1])
print("*** Plotting ODFs at strains:")
print(strain[ODF_tsteps])

#----------------------
# Initialize model
#----------------------

nlm_len  = sf.init(L)         # nlm_len is the number of fabric expansion coefficients (degrees of freedom).
lm       = sf.get_lm(nlm_len) # The (l,m) values corresponding to the coefficients in "nlm".
nlm      = np.zeros((Nt,nlm_len), dtype=np.complex128) # The expansion coefficients
nlm[0,0] = 1/np.sqrt(4*np.pi) # Init with isotropy (Normalized such that N(t=0) = 1)

#----------------------
# Integrate
#----------------------

vecdim = (Nt,3)
e1,e2,e3 = np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64),np.zeros(vecdim, dtype=np.float64)
eig = np.zeros(vecdim, dtype=np.float64)
Eij = np.zeros((Nt,3,3), dtype=np.float64)

# G=Glen, R=Rathmann & Lilien, P=Pettit, M=Martin
Y_G = sf.eps_of_tau__isotropic(tau0, Aglen,nglen)
Y_R, Y_P, Y_M = np.zeros((Nt,3,3), dtype=np.float64), np.zeros((Nt,3,3), dtype=np.float64), np.zeros((Nt,3,3), dtype=np.float64)

# Euler integration 
with Bar('dt=%.3fyr, Nt=%i :: L=%i (nlm_len=%i) ::'%(dt/year2sec,Nt,L,nlm_len), max=Nt-1, fill='#', suffix='%(percent).1f%% - %(eta)ds') as bar:
    for tt in np.arange(1,Nt+1):

        ttp = tt-1 # tt prev

        ### Prev state
        c  = nlm[ttp,:]

        ### Eigen frames (e_i)
        e1[ttp,:], e2[ttp,:], e3[ttp,:], eig[ttp,:] = sf.frame(c, 'e')

        ### Enhancement factors
        Eij[ttp,:,:] = np.transpose(sf.enhfac_eiej(c, e1[ttp,:],e2[ttp,:],e3[ttp,:], Ecc_lin, Eca_lin, alpha_lin, nprime))
        
        ### Y for fabric at constant strain rate
        Y_R[ttp,:,:] = sf.eps_of_tau__orthotropic(       tau0, Aglen,nglen, e1[ttp,:],e2[ttp,:],e3[ttp,:], Eij[ttp,:,:])
        Y_P[ttp,:,:] = sf.eps_of_tau__orthotropic_pettit(tau0, Aglen,nglen, e1[ttp,:],e2[ttp,:],e3[ttp,:], Eij[ttp,:,:])
        Y_M[ttp,:,:] = sf.eps_of_tau__orthotropic_martin(tau0, Aglen,nglen, e1[ttp,:],e2[ttp,:],e3[ttp,:], Eij[ttp,:,:])
 
        # Verify enhancement factors are correct
        if False:
            ei = np.zeros((3,3))
            ei[0,:] = e1[ttp,:]; ei[1,:] = e2[ttp,:]; ei[2,:] = e3[ttp,:];
            print("\n")
            print('e1,e2,e3 = ', ei[0,:],ei[1,:],ei[2,:])
            print("Eij = \n", Eij[ttp,:,:])
            for ii,jj in ((0,0),(1,1),(2,2), (1,2),(2,0),(0,1)):
                eij = np.tensordot(ei[ii,:],ei[jj,:],axes=0)
                tau0 = eij+eij.T if ii != jj else np.eye(3)/3-eij
                Eij_ = Eij[ttp,ii,jj]
                Y_G_ = np.tensordot(sf.eps_of_tau__isotropic(tau0, Aglen,nglen), eij, axes=2)
                Eij_R = np.tensordot(sf.eps_of_tau__orthotropic(       tau0, Aglen,nglen, e1[ttp,:],e2[ttp,:],e3[ttp,:], Eij[ttp,:,:]), eij, axes=2) / Y_G_
                Eij_P = np.tensordot(sf.eps_of_tau__orthotropic_pettit(tau0, Aglen,nglen, e1[ttp,:],e2[ttp,:],e3[ttp,:], Eij[ttp,:,:]), eij, axes=2) / Y_G_
                Eij_M = np.tensordot(sf.eps_of_tau__orthotropic_martin(tau0, Aglen,nglen, e1[ttp,:],e2[ttp,:],e3[ttp,:], Eij[ttp,:,:]), eij, axes=2) / Y_G_
                print("(i,j)=(%i,%i) :: Eij_R/Eij = %.3e --- Eij_P/Eij = %.3e --- Eij_M/Eij = %.3e :: Eij = %.3e"%(ii,jj, Eij_R/Eij_, Eij_P/Eij_, Eij_M/Eij_, Eij_))
            print("\n")
 
        if tt == Nt: break # We're done after calculating Y for the last step!
        
        ### Fabric evolution
        nu = sf.nu(nu0, eps)
        dndt = 0 * nu*sf.dndt_reg(c)
        
        if np.abs(strain[tt]) <= np.abs(DDRX_strainthres):
            dndt += nu*sf.dndt_reg(c)
            dndt += sf.dndt_latrot(c, eps, omg) 
        else:
            dndt += Gamma0 * sf.dndt_ddrx(c, tau0)  # TODO : should be tau0

        nlm[tt,:] = c + dt*np.matmul(dndt, c)
        
        bar.next()

#----------------------
# Plot
#----------------------

legkwargs = {'frameon':True, 'fancybox':False, 'edgecolor':'k', 'framealpha':1, 'ncol':1, 'handlelength':1.34, 'labelspacing':0.25, 'fontsize':FSSMALL}

color_darkred = '#a50f15'
color_blue = '#08519c'

color_R = 'k' 
color_P = 'k' 
color_M = 'k' 

### Setup figure

scale = 0.87

#fig = plt.figure(figsize=(8.0*scale,7.5*scale), constrained_layout=False)
#gs = gridspec.GridSpec(3, 2, height_ratios=[1,0.55,0.85])
#gs.update(hspace=0.44, wspace=0.26, left=0.1, right=0.96, top=0.98, bottom=0.08)

fig = plt.figure(figsize=(7.0*scale,5*scale), constrained_layout=False)
gs = gridspec.GridSpec(2, 2, height_ratios=[1,0.55])
gs.update(hspace=0.53, wspace=0.26, left=0.10, right=0.97, top=0.98, bottom=0.12)

ax_Y   = fig.add_subplot(gs[0, :])
#ax_E   = fig.add_subplot(gs[-1, 0])
#ax_eig = fig.add_subplot(gs[-1, 1])

gs_parcels = gridspec.GridSpecFromSubplotSpec(1, len(ODF_tsteps), subplot_spec=gs[1,:], width_ratios=[1,1,1,1], hspace=0.6, wspace=+0.2)
inclination = 50 # view angle
rot0 = -90
rot = -40 + rot0 
prj = ccrs.Orthographic(rot, 90-inclination)
geo = ccrs.Geodetic()     
ii0, jj0 = 0, 0
ax_ODF = [fig.add_subplot(gs_parcels[ii0, jj0+ii], projection=prj) for ii,nn in enumerate(ODF_tsteps)]

### Plot relative stresses

rect1 = plt.Rectangle((strain[0],0), DDRX_strainthres,  10, color='#deebf7')
rect2 = plt.Rectangle((DDRX_strainthres,0), strain[-1], 10, color='#fee0d2')
ax_Y.add_patch(rect1)
ax_Y.add_patch(rect2)

X = strain
if T_EXP==T_EXP_SS: xlims = X[[0,-1]] 
if T_EXP==T_EXP_CC: xlims = [0,-1]
if T_EXP==T_EXP_SS: strainMinorTicks = np.arange(0,20,0.5) 
if T_EXP==T_EXP_CC: strainMinorTicks = np.arange(-1,0,0.1)

if T_EXP==T_EXP_SS: 
    ii,jj = 0,2 # x,z
    ylbl = r'$\dot{\epsilon}_{xz}/\dot{\epsilon}^{\mathrm{Glen}}_{xz}$'
    xglen,yglen = 4.5,1.05
    
if T_EXP==T_EXP_CC: 
    ii,jj = 2,2 # z,z
    ylbl = r'$\dot{\epsilon}_{zz}/\dot{\epsilon}^{\mathrm{Glen}}_{zz}$'
    xglen,yglen = -0.625,1.08

cglen = '0.4'
ax_Y.plot(X, X*0+1, '-', lw=1.2,  color=cglen)
ax_Y.text(xglen, yglen, "Glen's law", color=cglen, horizontalalignment='center', fontsize=FSSMALL)

Yn_R = Y_R[:,ii,jj]/Y_G[ii,jj]
Yn_P = Y_P[:,ii,jj]/Y_G[ii,jj]
Yn_M = Y_M[:,ii,jj]/Y_G[ii,jj]
delta_Martin = 100*np.divide(np.abs(Yn_M-Yn_R),Yn_R)
delta_Pettit = 100*np.divide(np.abs(Yn_P-Yn_R),Yn_R)
print('max(Martin-True)=%.2f'%(np.amax(delta_Martin)))
print('max(Pettit-True)=%.2f'%(np.amax(delta_Pettit)))
print('(Martin-True)_{-1}=%.2f'%(delta_Martin[-1]))
print('(Pettit-True)_{-1}=%.2f'%(delta_Pettit[-1]))
Yall = np.vstack((Yn_R,Yn_P,Yn_M))
ax_Y.plot(X, Yn_R, '-',  color=color_R, label=r"Unapprox.")
ax_Y.plot(X, Yn_M, '--', color=color_M, label=r"Martin")
ax_Y.plot(X, Yn_P, ':',  color=color_P, label=r"Pettit")
ax_Y.set_xlabel(r'$\epsilon$')
ax_Y.set_xticks(strainMinorTicks, minor=True)  
ax_Y.set_xlim(xlims)
ax_Y.set_ylabel(ylbl)

leg = ax_Y.legend(loc=4, **legkwargs) # bbox_to_anchor=(1,0.88)
leg.get_frame().set_linewidth(0.8);
writeSubplotLabel(ax_Y,2,r'{\bf a}')

dylim = 0.07*(np.amax(Yall)-np.amin(Yall)) # % of data range
ylims = [np.amin(Yall)-dylim, np.amax(Yall)+dylim]

if T_EXP==T_EXP_SS:
    ax_Y.set_yticks(np.arange(0,10+1e-3,0.5))
    ax_Y.set_yticks(np.arange(0,10+1e-3,0.1), minor=True)
    
if T_EXP==T_EXP_CC:
    ax_Y.set_yticks(np.arange(0,10+1e-3,0.5))
    ax_Y.set_yticks(np.arange(0,10+1e-3,0.1), minor=True)

ax_Y.set_ylim(ylims)
dx, y0 = 0.016, ylims[-1]-1.1*dylim
ax_Y.text(DDRX_strainthres*(1-dx), y0, r'{\bf Lattice rotation}', color=color_blue,    horizontalalignment='right', verticalalignment='center', fontsize=FSSMALL)
ax_Y.text(DDRX_strainthres*(1+dx), y0, r'{\bf DDRX}',             color=color_darkred, horizontalalignment='left', verticalalignment='center', fontsize=FSSMALL)

#### Plot enahancement factors          

#ax_E.plot(X, Eij[:,0,2],  '--', color='k',  label=r'$E_{13}$')
#ax_E.plot(X, Eij[:,0,1],  ':',  color='k',  label=r'$E_{12}$')
#ax_E.plot(X, Eij[:,0,0],  '-',  color='k',  label=r'$E_{11}$')
#ax_E.set_xlabel(r'$\epsilon$')
#ax_E.set_xticks(strainMinorTicks, minor=True)  
#ax_E.set_xlim(xlims)
#ax_E.set_ylabel(r'$E_{ij}$')
#ylims = [0,2.5]
#ax_E.set_yticks(np.arange(ylims[0],ylims[-1]+1e-4,1))
#ax_E.set_yticks(np.arange(ylims[0],ylims[-1]+1e-4,0.5), minor=True)
#ax_E.set_ylim(ylims)
##ax_E.grid()
#leg = ax_E.legend(loc=3 if T_EXP==T_EXP_CC else 5, **legkwargs)
#leg.get_frame().set_linewidth(0.8);
#writeSubplotLabel(ax_E,2,r'{\bf f}')

#### Plot fabric eigenvalues

#ax_eig.plot(X, eig[:,0],  '-',  color='k',  label=r'$\lambda_1$')
#ax_eig.plot(X, eig[:,1],  '--', color='k',  label=r'$\lambda_2$')
#ax_eig.plot(X, eig[:,2],  ':',  color='k',  label=r'$\lambda_3$')
#ax_eig.set_xlabel(r'$\epsilon$')
#ax_eig.set_xticks(strainMinorTicks, minor=True)  
#ax_eig.set_xlim(xlims)
#ax_eig.set_ylabel(r'$\lambda_i$')
#ax_eig.set_yticks([0,0.5,1])
#ax_eig.set_yticks(np.arange(0,1.1,0.1), minor=True)
#ax_eig.set_ylim([0,1])
##ax_eig.grid()
#ax_eig.legend(loc=1 if T_EXP==T_EXP_CC else 5, **legkwargs)
#writeSubplotLabel(ax_eig,2,r'{\bf g}')

### Plot ODFs

colax = 'k'

def plot_ei(ax, vec, geo, ei_num, mrk='.', ms=7, lbl=False, phi_rel=0, theta_rel=0):
    theta,phi = getPolarAngles(vec)
    if 0<=phi<=180 and (theta<=0): theta,phi = getPolarAngles(-vec)
    color = color_darkred
    ax.plot([phi],[theta], mrk, ms=ms, markerfacecolor=color, markeredgecolor=color, markeredgewidth=1.0, transform=geo)

def getPolarAngles(vec):
    x,y,z = vec
    phi   = np.rad2deg(np.arctan2(y,x))
    theta = 90 - np.rad2deg(np.arccos(z))
    return (theta, phi)

for ii,nn in enumerate(ODF_tsteps):

    ax = ax_ODF[ii]
    ax.set_global()
    N = np.sqrt(4*np.pi)*nlm[nn,0]
    ODF = np.divide(nlm[nn,:],N) # normalize distribution (=ODF)
    plot_ODF(ODF, lm, ax=ax, cmap='Greys', cblabel=r'ODF')
    ax.set_title(r'$\epsilon=%.1f$'%(strain[nn]), fontsize=FS)
    writeSubplotLabel(ax,2,r'{\bf %s}'%(chr(ord('b')+ii)), bbox=(-0.25,1.35))

    if ii==0:
        mrk='X'
        ms=4.5
        ax.plot([0],[90],mrk, ms=ms, c=colax, transform=geo)
        ax.plot([rot0-90],[0],mrk, ms=ms,  c=colax, transform=geo)
        ax.plot([rot0],[0],mrk, ms=ms,  c=colax, transform=geo)
        ax.text(rot0-80, 65, r'$\vu{z}$', c=colax, horizontalalignment='left', transform=geo)
        ax.text(rot0-90+12, -5, r'$\vu{x}$', c=colax, horizontalalignment='left', transform=geo)
        ax.text(rot0-23, -5, r'$\vu{y}$', c=colax, horizontalalignment='left', transform=geo)
    
    if ii>0:    
        color, ms = color_darkred, 7
#        w,v = np.linalg.eig(tau_RL[nn,:,:]) # debug: principal stress directions
#        for ei in (v[:,0],v[:,1],v[:,2]): # debug: principal stress directions
        for ei in (e1[nn,:],e2[nn,:],e3[nn,:]):
            theta,phi = getPolarAngles(ei)
            ax.plot([phi],[theta], marker='.', ms=ms, markerfacecolor=color, markeredgecolor=color, markeredgewidth=1.0, transform=geo)
            theta,phi = getPolarAngles(-ei)
            ax.plot([phi],[theta], marker='.', ms=ms, markerfacecolor=color, markeredgecolor=color, markeredgewidth=1.0, transform=geo)

    if False and ii==0:    
        color, ms = color_blue, 8
        w,v = np.linalg.eig(tau0) # debug: principal stress directions
        #print(w,v)
        for ei in (v[:,0],v[:,1],v[:,2]): # debug: principal stress directions
            ei = np.array(ei)
            theta,phi = getPolarAngles(ei)
            ax.plot([phi],[theta], marker='o', ms=ms, markerfacecolor='none', markeredgecolor=color, markeredgewidth=1.0, transform=geo)
            theta,phi = getPolarAngles(-ei)
            ax.plot([phi],[theta], marker='o', ms=ms, markerfacecolor='none', markeredgecolor=color, markeredgewidth=1.0, transform=geo)

        
### Plot parcel deformations (on seperated "inset" axes)

# Parcel geometry
xyz0_init = (1,1,1)
ex,ey,ez = np.array([xyz0_init[0],0,0]),np.array([0,xyz0_init[1],0]),np.array([0,0,xyz0_init[2]])
dzx, dzy, dyx = 0,0,0
time_ = dt*np.arange(0,Nt) 

axsize = 0.14
scale = 1

if T_EXP==T_EXP_SS: axy0 = 0.70
if T_EXP==T_EXP_CC: axy0 = 0.70

if T_EXP==T_EXP_SS: steps_to_plot_ODF = ODF_tsteps[0:2]
if T_EXP==T_EXP_CC: steps_to_plot_ODF = [ODF_tsteps[0], ODF_tsteps[2]] 

for ii,nn in enumerate(steps_to_plot_ODF):

    if T_EXP==T_EXP_SS:
        dyx = np.dot(np.matmul(SS.F(time_[nn])-np.eye(3), ey), ex)
        dzx = np.dot(np.matmul(SS.F(time_[nn])-np.eye(3), ez), ex)
        dzy = 0
        pcoords = np.array([[0.106,axy0],[0.28,axy0]])
        
    if T_EXP==T_EXP_CC:
        xyz0_init = np.matmul(PS.F(time_[nn]), np.array(xyz0_init))
        pcoords = np.array([[0.106,axy0],[0.55,axy0]])

    ax_sub = fig.add_axes([pcoords[ii,0],pcoords[ii,1], axsize,axsize], projection='3d')
    ax_sub.patch.set_alpha(0.01)
    plot_parcel(ax_sub, xyz0_init, dzx,dzy,dyx, color='0.5', plotaxlbls=True, scale=scale)
    ax_sub.set_title(r'$\epsilon=%.1f$'%(strain[nn]), fontsize=FSSMALL,  x=0.5, y=0.92)
    
    if T_EXP==T_EXP_SS:
        
        lw0 = 1.25
        lwdotted = 0.65* lw0
        coldotted = '0.4'
        x0, *_ = xyz0_init
        zero, one, doubleone = np.array([0,0]), np.array([0,1*scale]), np.array([1*scale,1*scale])
        zspan = one
        ax_sub.plot(np.array([0,dzx])*scale,zero,doubleone, ':', lw=lwdotted, color=coldotted, zorder=10, clip_on=False)
        ax_sub.plot(np.array([x0+dzx]*2)*scale,zero,zspan,  ':', lw=lwdotted, color=coldotted, zorder=10, clip_on=False)
        ax_sub.plot(np.array([x0,x0+dzx])*scale,zero,zero,  ':', lw=lwdotted, color=coldotted, zorder=10, clip_on=False)

### Save figure

fname = 'rheology-compare-%s.png'%(T_EXP_STR[T_EXP])
print('Saving output to %s'%(fname))
plt.savefig(fname, dpi=300)

