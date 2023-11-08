# N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

import sys, os, copy, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import plotting as sfplt
FS = sfplt.setfont_tex(fontsize=13)

import warnings
warnings.filterwarnings("ignore")

### Run options

MAKE_FRAMES = True
MAKE_ANI    = 0

Mtypes = ['LROT','DDRX','CDRX']
#Mtypes = ['LROT','CDRX']

### Numerics

Nt = 100 # Number of time steps
L  = 8

### Integrate

lm, nlm_len = sf.init(L) 
nlm = {\
    'LROT': np.zeros((3,Nt+1,nlm_len), dtype=np.complex64), \
    'DDRX': np.zeros((3,Nt+1,nlm_len), dtype=np.complex64), \
    'CDRX': np.zeros((3,Nt+1,nlm_len), dtype=np.complex64), \
}

for ii, Mtype in enumerate(Mtypes):
    for jj, exptype in enumerate(['uc_zz', 'cc_zy', 'ss_yz']):

        # Process rate magnitude
        if Mtype == 'LROT': kwargs = dict(iota=1,    Lambda=None, Gamma0=None, nu=1)
        if Mtype == 'CDRX': kwargs = dict(iota=None, Lambda=2e-1, Gamma0=None, nu=1)
        if Mtype == 'DDRX': kwargs = dict(iota=None, Lambda=None, Gamma0=2e1,  nu=1)

        # Initial state
        if Mtype == 'CDRX': nlm0 = nlm['LROT'][jj,-1,:]
        else:               nlm0 = None
        
        # Mode of deformation
        if exptype == 'uc_zz': # Uniaxial compression
            mod = dict(type='ps', r=0,  axis=2)
            strain_target = -0.90

        if exptype == 'cc_zy': # Confined compression
            mod = dict(type='ps', r=-1, axis=2)
            strain_target = -0.90

        if exptype == 'ss_yz': # Simple shear
            mod = dict(type='ss', plane=0)
            strain_target = np.deg2rad(75) if Mtype == 'LROT' else np.deg2rad(66.5) # => dt as uc/cc if not LROT
        
        # Integrate 
        nlm[Mtype][jj,:,:], *_ = sfint.lagrangianparcel(sf, mod, strain_target, Nt=Nt, nlm0=nlm0, **kwargs)

#Mtypes = ['CDRX',]
print('Model integrations finished...')
        
### Plot

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

geo, prj = sfplt.getprojection(rotation=55-90, inclination=50)

for ii, Mtype in enumerate(Mtypes):

    if MAKE_FRAMES: 
    
        for tt in np.arange(0,Nt):
        #for tt in np.arange(Nt-1,Nt): # debug
        
            ### Figure setup
            
            scale = 2.6
            fig = plt.figure(figsize=(3/2*1.6*scale,1.00*scale))
            gs = fig.add_gridspec(1,3)
            gs.update(left=0.04, right=1-0.02, top=0.85, bottom=0.20, wspace=0.3, hspace=0.4)

            axi = [fig.add_subplot(gs[0,jj], projection=prj) for jj in range(3)]
            for jj in range(3): axi[jj].set_global(); 

            ### Plot
            
            lvlset = [np.linspace(0,1,9), lambda x,p:'%.1f'%x]
            for jj in range(3): sfplt.plotODF(nlm[Mtype][jj,tt,:], lm, axi[jj], lvlset=lvlset)
            for jj in range(3): sfplt.plotcoordaxes(axi[jj], geo, axislabels='vuxi')
            
            axi[0].set_title(r'%s Unconfined pure shear'%(r'{\Large \textit{(a)}}\,\,'), fontsize=FS, pad=10)
            axi[1].set_title(r'%s Confined pure shear'%(r'{\Large \textit{(b)}}\,\,'), fontsize=FS, pad=10)
            axi[2].set_title(r'%s Simple shear'%(r'{\Large \textit{(c)}}\,\,'), fontsize=FS, pad=10)
                
            ### Save
            
            fout = 'frames/fabric-evolution-%s-animation-%03i.png'%(Mtype,tt)
            print('Saving %s'%(fout))
            plt.savefig(fout, dpi=200)
            plt.close()

    if MAKE_ANI:
    
        os.system('rm animation-%s.gif'%(Mtype))  
        os.system('ffmpeg -y -f image2 -framerate 50 -stream_loop 0 -i frames/fabric-evolution-%s-animation-%s.png -vcodec libx264 -crf 20  -pix_fmt yuv420p animation-%s.avi'%(Mtype,'%03d',Mtype))
        os.system('ffmpeg -i animation-%s.avi -vf "fps=20,scale=550:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 animation-%s.gif'%(Mtype,Mtype))
        os.system('rm animation-%s.avi'%(Mtype))  

