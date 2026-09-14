# N. M. Rathmann <rathmann@nbi.ku.dk>, 2026-

import sys, os, copy, code # code.interact(local=locals())
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
from matplotlib import rcParams, rc
import matplotlib.transforms as mtransforms

import cmasher as cmr

sys.path.insert(0, '..')
#from localheader import *

from specfabpy import specfab as sf
from specfabpy import common as sfcom
from specfabpy import discrete as sfdsc
from specfabpy import plotting as sfplt
from specfabpy import integrator as sfint

### Fonts

FS = 11
FS = sfplt.setfont_tex(fontsize=FS)
FSSMALL = FS-1

rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage[sc]{mathpazo} \usepackage{siunitx} \DeclareSIUnit\year{a}'


### Construct states

L = 14
lm, nlm_len = sf.init(L) 

# Ideal

#nlm_iso    = sfcom.nlm_isotropic(L)    # isotropic
nlm_pole   = sf.nlm_idealz(0, L)       # delta pole
nlm_girdle = sf.nlm_idealz(np.pi/2, L) # delta girdle
#nlm_circle = sf.nlm_idealz(np.deg2rad(70), L) # delta circle

Sl_pole,   Lrng = sf.powerspectrum(nlm_pole, L)
Sl_girdle, Lrng = sf.powerspectrum(nlm_girdle, L)
#Sl_circle, Lrng = sf.powerspectrum(nlm_circle, L)

# Weaker versions

kwargs = dict(Nt=100, iota=+1)

strain_target = -0.8
nlm_wpole,   F, *_ = sfint.lagrangianparcel(sf, dict(type='ps', axis=2, tau=+1, q=0) , strain_target, **kwargs, nu=1)

strain_target = 6
nlm_wgirdle, F, *_ = sfint.lagrangianparcel(sf, dict(type='ps', axis=2, tau=-1, q=0) , strain_target, **kwargs, nu=1)  

nlm_wpole   = nlm_wpole[-1]
nlm_wgirdle = nlm_wgirdle[-1]

Sl_wpole, _   = sf.powerspectrum(nlm_wpole, L)
Sl_wgirdle, _ = sf.powerspectrum(nlm_wgirdle, L)

### Plot

scale=0.6
fig = plt.figure(figsize=(4.7*scale,3.5*scale))
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0,0])

lw = 1.25
kw = dict(lw=lw)
kww = dict(lw=lw)
cp, cg = "#67000d", "#000000"
cwp, cwg = "#cb181d", "#737373"
#cp, cg = "#000000", "#000000"
#cwp, cwg = "#000000", "#000000"

hp,  = ax.semilogy(Lrng, Sl_pole,    c=cp,  label=r'Perfect pole', ls='-', **kw)
hwp, = ax.semilogy(Lrng, Sl_wpole,   c=cwp, label=r'Weaker pole', ls=':', **kww)
hg,  = ax.semilogy(Lrng, Sl_girdle,  c=cg,  label=r'Perfect girdle', ls='--', **kw)
hwg, = ax.semilogy(Lrng, Sl_wgirdle, c=cwg, label=r'Weaker girdle', ls='-.', **kww)

# Set figure axes etc.

ax.set_xlabel('$l$')
ax.set_ylabel('$S_l/S_0$')
ax.set_ylim([4e-3,1.5])
ax.set_xticks(np.arange(0,21,2))
#ax.set_xticks(np.arange(0,21,2), minor=True)
ax.set_xlim([0,np.amax([10,L])])
#ax.grid()

#legkwargs = {'frameon':False, 'fancybox':False, 'edgecolor':'k', 'columnspacing':1, 'handlelength':1.2, 'labelspacing':0.2}
##hleg = ax.legend(loc=3, fontsize=FS, ncol=2, bbox_to_anchor=(0.0,-0.45), **legkwargs)
#hleg = ax.legend(loc=3, fontsize=FS-0.5, ncol=1, bbox_to_anchor=(-0.01,-0.04), **legkwargs)

def label_line(ax, x, y, label, x_label, offset_points=-3, window=2, ha='left', color=None, **text_kwargs):

    x, y  = np.asarray(x), np.asarray(y)

    i = int(np.argmin(np.abs(x - x_label))) # index closest to the requested x position

    # clamp the window so we don't run off the ends of the array
    lo = max(i - window, 0)
    hi = min(i + window, len(x) - 1)
    if lo == hi:  # degenerate: only one point available
        lo, hi = max(i - 1, 0), min(i + 1, len(x) - 1)

    # slope in display (pixel) space so log axes and aspect are handled
    p1 = ax.transData.transform((x[lo], y[lo]))
    p2 = ax.transData.transform((x[hi], y[hi]))
    dx_disp, dy_disp = p2 - p1
    angle = np.degrees(np.arctan2(dy_disp, dx_disp))

    # offset in points, so the gap is constant regardless of scale/zoom
    offset = mtransforms.offset_copy(ax.transData, fig=ax.figure, x=-3, y=offset_points, units='points')

    va = 'top' if offset_points < 0 else 'bottom'

    return ax.text(x[i], y[i], label, transform=offset, ha=ha, va=va, rotation=angle, rotation_mode='anchor', color=color, fontsize=FSSMALL, **text_kwargs)
    
                   
label_line(ax, Lrng, Sl_pole,    "Perfect pole",   x_label=14,  offset_points=-3, color=hp.get_color(), ha='right')
label_line(ax, Lrng, Sl_wpole,   "Weaker pole",    x_label=5,  offset_points=+2, color=hwp.get_color())
label_line(ax, Lrng, Sl_girdle,  "Perfect girdle", x_label=14,  offset_points=+3, color=hg.get_color(), ha='right')
label_line(ax, Lrng, Sl_wgirdle, "Weaker girdle",  x_label=6,  offset_points=-6, color=hwg.get_color(), ha='center')

#i = len(Lrng) // 1.2
#offset = mtransforms.offset_copy(ax.transData, fig=ax.figure, x=0, y=-3, units='points')

#ax.text(Lrng[i], Sl_pole[i], "Perfect pole", transform=offset, ha='center', va='top')

#print('pfJ (new / old) for iso',    sf.pfJ(nlm_iso, L), sfcom.pfJ(nlm_iso))
print('pfJ (new / old) for pole',   sf.pfJ(nlm_pole, L), sfcom.pfJ(nlm_pole))
print('pfJ (new / old) for girdle', sf.pfJ(nlm_girdle, L), sfcom.pfJ(nlm_girdle))

### Plot ODF insets?

if 1:

    geo, prj = sfplt.getprojection(rotation=55+180, inclination=50)

    x0 = 0.92
    y0 = 0.02
    dy = 0.22
    ODF_plots = (\
        {'nlm':nlm_pole,    'title':'Perfect pole',   'axloc':(x0, y0+3*dy), 'lvlmax':1.0, 'lvlmin':0.6}, \
        {'nlm':nlm_wpole,   'title':'Weaker pole',    'axloc':(x0, y0+2*dy), 'lvlmax':0.4, 'lvlmin':0.15}, \
        {'nlm':nlm_girdle,  'title':'Perfect girdle', 'axloc':(x0, y0+1*dy), 'lvlmax':0.6, 'lvlmin':0.15}, \
        {'nlm':nlm_wgirdle, 'title':'Weaker girdle',  'axloc':(x0, y0+0*dy), 'lvlmax':0.7, 'lvlmin':0.15}, \
    )

    for ODF in ODF_plots:

        W = 0.19 # ax width
        axpos = [ODF['axloc'][0], ODF['axloc'][1], W,W]
        axin = plt.axes(axpos, projection=prj) #, transform=ax.transData)
        axin.set_global()
        
        nlm = ODF['nlm']
        cmap = cmr.get_sub_cmap('Reds' if "pole" in ODF['title'] else 'Greys', 0.25, 1) # don't include pure white.
        cmap.set_under('w')
        lvlset = [np.linspace(ODF['lvlmin'], ODF['lvlmax'], 8), lambda x,p:'%.1f'%x]
        sfplt.plotODF(nlm, lm, axin, lvlset=lvlset, cmap=cmap, showcb=False)

        axin.text(1.15, 0.5, ODF['title'], transform=axin.transAxes, fontsize=FSSMALL, ha='left', va='center', clip_on=False) 

### Save figure

plt.savefig('power-spectrum.png', bbox_inches='tight', pad_inches=0.02, dpi=300)

