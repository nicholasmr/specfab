# Lattice rotation

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/plastic-spin.png){: style="width:350px"}

When subject to simple shear, the orientation of easy slip systems in polycrystalline materials like ice or olivine tend towards aligning with the bulk shear-plane system. 
That is, the $n$-axis ($c$-axis in ice) tends towards to bulk shear plane normal, and the $b$-axis ($a$-axis in ice) towards to bulk shear direction.
Thus, if grains are perfectly aligned with the bulk shear system, their orientation should be unaffected by any further shearing, but be in steady state.
Clearly, slip systems do therefore not simply co-rotate with the bulk continuum spin ($\bf W$) like passive material line elements embedded in a flow field, i.e. 
${\bf \dot{n}} \neq {\bf W} \cdot {\bf n}$.
Rather, slip systems must be subject to an additional contribution &mdash; a plastic spin ${\bf W}_\mathrm{p}$ &mdash; such that the bulk spin is exactly counteracted to achieve steady state if favorably aligned:

$$
{\bf \dot{n}} = ({\bf W} + {\bf W}_{\mathrm{p}}^{(n)}) \cdot {\bf n} = {\bf 0} \quad\text{for ${\bf b}$&ndash;${\bf n}$ shear}.
$$

$$
{\bf \dot{b}} = ({\bf W} + {\bf W}_{\mathrm{p}}^{(b)}) \cdot {\bf b} = {\bf 0} \quad\text{for ${\bf b}$&ndash;${\bf n}$ shear}.
$$

More precisely, the crystallographic axes reorient themselves in response to both the bulk continuum spin and a plastic spin that is supposed to represent the crystallographic spin needed to accommodate strain compatibility between grains that otherwise preferentially deform by easy slip (basal slip for ice).

## Directors method

[Aravas and Aifantis (1991)](https://doi.org/10.1016/0749-6419(91)90028-W) and [Aravas (1994)](https://www.doi.org/10.1088/0965-0393/2/3A/005) (among others) proposed a particularly simple model for the functional form of ${\bf W}_\mathrm{p}$, the so-called directors method.

For a constant rate of shear deformation ($1/T$) aligned with the ${\bf b}$&mdash;${\bf n}$ system, 

$$
{\bf F} = {\bf I} + \frac{t}{T} {\bf b}\otimes{\bf n}
\quad \Longrightarrow\quad
\nabla {\bf u} = \frac{1}{T} {\bf b}\otimes{\bf n}
,
$$

the bulk strain-rate and spin tensors are, respectively,

$$
{\bf D} = \frac{1}{2T} ({\bf b}\otimes{\bf n} + {\bf n}\otimes{\bf b}),
\\
{\bf W} = \frac{1}{2T} ({\bf b}\otimes{\bf n} - {\bf n}\otimes{\bf b}) .
$$

Since ${\bf W}_\mathrm{p} = -{\bf W}$ is required in steady state, it follows from eliminating $1/(2T)$ by calculating 
${\bf D} \cdot {\bf n}\otimes{\bf n}$, 
${\bf n}\otimes{\bf n} \cdot {\bf D}$, 
${\bf D} \cdot {\bf b}\otimes{\bf b}$, and 
${\bf b}\otimes{\bf b} \cdot {\bf D}$, 
that

$$
{\bf W}_\mathrm{p}^{(n)} = +{\bf n}\otimes{\bf n} \cdot {\bf D} - {\bf D} \cdot {\bf n}\otimes{\bf n} ,
\\
{\bf W}_\mathrm{p}^{(b)} = -{\bf b}\otimes{\bf b} \cdot {\bf D} + {\bf D} \cdot {\bf b}\otimes{\bf b} ,
$$

Indeed, this result agrees with representation theorems for isotropic functions (Wang, 1969), stating that an antisymmetric tensor-valued function of a symmetric tensor (${\bf D}$) and a vector (${\hat {\bf r}}$) is to lowest order given by

$$
{\bf W}_{\mathrm{p}} = 
\iota({\hat {\bf r}}\otimes{\hat {\bf r}}\cdot{\bf D} - {\bf D}\cdot{\hat {\bf r}}\otimes{\hat {\bf r}})
.
$$

To be consistent with the above, $\iota = +1$ for ${\hat {\bf r}}={\bf n}$ and $\iota = -1$ for ${\hat {\bf r}}={\bf b}$.

${\bf W}_\mathrm{p}^{(n)}$ and ${\bf W}_\mathrm{p}^{(b)}$ are then generally taken to suffice for other deformation kinematics, too.

!!! note "Glacier ice"

    The predicted normalized $c$-axis velocity field for glacier ice (${\bf \dot{c}}={\bf \dot{n}}$) is show below for three different deformation kinematics:

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation-modes/deformation-modes.png#center){: style="width:550px"}

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/latrot-velocity/latrot-velocity.png){: style="width:560px"}

### Matrix model 

The corresponding effect on the continuous distribution functions is modelled as a conservative advection process in orientation space $S^2$:

$$ 
\frac{\mathrm{D} n}{\mathrm{D} t} = -\nabla_{S^2}\cdot(n{\bf \dot{n}}) 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}_n}{\mathrm{D} t} = {\bf M_{\mathrm{LROT}}}(\iota=+1) \cdot {\bf s}_n,
$$

$$ 
\frac{\mathrm{D} b}{\mathrm{D} t} = -\nabla_{S^2}\cdot(b{\bf \dot{b}}) 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}_b}{\mathrm{D} t} = {\bf M_{\mathrm{LROT}}}(\iota=-1) \cdot {\bf s}_b,
$$

where ${\bf M_{\mathrm{LROT}}}(\iota)$ is given analytically in [Rathmann et al. (2021)](https://doi.org/10.1017/jog.2020.117).

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/fabric-evolution/animation-LROT.gif){: style="width:650px"}


### Code example

```python
import numpy as np
from specfabpy import specfab as sf
# L=8 truncation is sufficient in this case, but larger L allows a very strong fabric to  
#  develop and minimizes the effect that regularization has on low wavenumber modes (l=2,4)
lm, nlm_len = sf.init(8) 

### Velocity gradient tensor experienced by parcel

ugrad = np.diag([0.5, 0.5, -1.0]) # uniaxial compression along z-axis
D = (ugrad+np.transpose(ugrad))/2 # symmetric part (strain rate tensor)
W = (ugrad-np.transpose(ugrad))/2 # anti-symmetric part (spin tensor)

### Numerics 

Nt = 25   # number of time steps
dt = 0.05 # time-step size

### Initial fabric state

nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # n state vector
nlm[0,0] = 1/np.sqrt(4*np.pi) # normalized isotropic state at t=0

### Euler integration

for tt in np.arange(1,Nt):
    nlm_prev = nlm[tt-1,:] # previous solution
    iota, zeta = 1, 0      # "deck of cards" behavior 
    M_LROT = sf.M_LROT(nlm_prev, D, W, iota, zeta) # lattice rotation operator
    M_REG  = sf.M_REG(nlm_prev, D)                 # regularization operator
    M      = M_LROT + M_REG
    nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # euler step
    nlm[tt,:] = sf.apply_bounds(nlm[tt,:]) # apply spectral bounds if needed
```

