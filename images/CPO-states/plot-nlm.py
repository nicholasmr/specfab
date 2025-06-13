#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2023-

import numpy as np
import copy, sys, os, code # code.interact(local=locals())
import cartopy.crs as ccrs

from specfabpy import specfab as sf
from specfabpy import plotting as sfplt

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

lm, nlm_len = sf.init(4) 
FS = sfplt.setfont_tex(fontsize=13)
geo, prj = sfplt.getprojection(rotation=-90-45, inclination=45)
kw_fig = dict(transparent=True, bbox_inches='tight', pad_inches=0.025)
 
def plot_nlm(a2, fname, ei=None, title='n', micolor='black'):

    nlm = np.zeros(nlm_len, dtype=np.complex128)
    nlm[:sf.L2len] = sf.a2_to_nlm(a2)

    scale = 1
    fig = plt.figure(figsize=(1.125*scale,1.3*scale))
    gs = gridspec.GridSpec(1, 1)
    
    a=0.14
    gs.update(top=0.80, bottom=0.00, left=a, right=1-a)
    ax = fig.add_subplot(gs[0, 0], projection=prj)
    ax.set_global()
    if title is not None: ax.set_title(r'$%s(\theta,\phi)$'%(title), pad=7, fontsize=FS)

    lvlset = (np.linspace(0,0.30,9), lambda x,p:'%.1f'%x) 
    nchunk = None if np.sum(nlm[1:])> 1e-2 else 49
    sfplt.plotODF(nlm, lm, ax, lvlset=lvlset, nchunk=nchunk, showcb=False)
    
    kw_axes = dict(transform=geo, color=micolor, fontsize=FS)
    if ei=='mt':
        sfplt.plotS2text(ax, [0,0,1],  r'$\vb{m}$', **kw_axes)
        sfplt.plotS2text(ax, [0,-1,0], r'$\vb{t}$', **kw_axes)
    elif ei=='mi': 
        sfplt.plotS2text(ax, [0,0,1],  r'$\vb{m}_1$', **kw_axes)
        sfplt.plotS2text(ax, [0,-1,0], r'$\vb{m}_2$', **kw_axes)
        sfplt.plotS2text(ax, [-1,0,0], r'$\vb{m}_3$', **kw_axes)    
    elif ei=='mialt': 
        # radioglacology convention where \lambda_3 is the largest eigenvalue (as opposed to \lambda_1)
        sfplt.plotS2text(ax, [0,0,1],  r'$\vb{m}_3$', **kw_axes)
        sfplt.plotS2text(ax, [0,-1,0], r'$\vb{m}_2$', **kw_axes)
        sfplt.plotS2text(ax, [-1,0,0], r'$\vb{m}_1$', **kw_axes)    
        
    fields = [fname,]
    if ei is not None:     fields.append(ei)
    if micolor != 'black': fields.append('red')
    if title is None:      fields.append('notitle')
    
    plt.savefig('nlm-%s.pdf'%('-'.join(fields)), **kw_fig)
    plt.close()
    return ax
    
### States 

states = dict(
    iso      = {'a2':np.eye(3)/3,              'ei':(None,'mi','mialt')},
    traniso  = {'a2':np.diag([0.2, 0.2, 0.6]), 'ei':(None,'mt')}, 
    ortho    = {'a2':np.diag([0.1, 0.3, 0.6]), 'ei':(None,'mi')}, 
    smaxz    = {'a2':np.diag([0.2, 0.2, 0.6]), 'ei':(None,'mi','mialt')}, 
    girdlexz = {'a2':np.diag([0.2, 0.4, 0.4]), 'ei':(None,'mi','mialt')}, 
)

for name, state in states.items():
    for ei in state['ei']:
        for c in ('black', sfplt.c_dred):
            if ei is None and c != 'black': continue
            plot_nlm(state['a2'], name, ei=ei, micolor=c)
            plot_nlm(state['a2'], name, ei=ei, micolor=c, title=None)
     
