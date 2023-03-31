# CPO Dynamics

## Introduction

This tutorial shows how to model the CPO evolution of a Lagrangian material parcel.<br>
Three modes of deformation/stress are considered:

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation-modes/deformation-modes.png#center){: style="width:620px"}

### Notation

Given the expansion

$$
n(\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}n_{l}^{m}Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip-plane normals)}, 
$$

CPO evolution can be written as a matrix problem involving

$$
{\boldsymbol n} = [n_0^0,n_2^{-2},n_2^{-1},n_2^{0},n_2^{1},n_2^{2},n_4^{-4},\cdots,n_4^{4},\cdots,n_L^{-L},\cdots,n_L^{L}]^{\mathrm{T}} \quad\text{(state vector)},
$$

such that 

$$
\frac{\partial}{\partial t} {\boldsymbol n} = {\bf M} \cdot {\boldsymbol n} \quad\text{(state evolution)},
$$

where the operator (matrix) ${\bf M}$ represents the effect of a given CPO process, which may depend on stress, strain-rate, temperature, etc.

The total effect of multiple processes acting simultaneously is simply

$$
{\bf M} = {\bf M_{\mathrm{LROT}}} + {\bf M_{\mathrm{DDRX}}} + {\bf M_{\mathrm{CDRX}}} + \cdots \quad\text{(operator)}. 
$$

!!! note
    The tutorial focuses only on modelling the CPO evolution of glacier ice, i.e. the $c$-axis (or mass fraction) distribution $n(\theta,\phi)$.

    $n(\theta,\phi)$ is also referred to as $\psi(\theta,\phi)$ in some of the figures below. 

- - -

## Lattice rotation

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/fabric-evolution/animation-LROT.gif){: style="width:650px"}

The strain-induced rotation of $c$-axes is modelled given the (local) velocity gradient tensor $\nabla {\bf u}$ ([Svendsen and Hutter, 1996](https://doi.org/10.3189/S0260305500013525)).
The model is a kinematic model in the sense that $c$-axes rotate in response to the bulk rate of stretching, ${\bf D}$, and spin, ${\bf W}$, thereby allowing the detailed microscopic stress and strain rate fields to be neglected and hence interactions between neighboring grains to be disregarded.

The modelled $c$-axes are taken to rotate with the bulk continuum spin (${\bf W}$), plus some plastic spin correction (${\bf W}_{\mathrm{p}}$), so that the $c$-axis velocity field on the unit sphere is

$$
{\bf \dot{c}} = ({\bf W} + {\bf W}_{\mathrm{p}}) \cdot {\bf c} \quad\text{($c$-axis velocity field)}
,
$$

where the plastic spin generally depends on ${\bf D}$ to lowest and second lowest order as (Wang, 1969; [Aravas, 1994](https://www.doi.org/10.1088/0965-0393/2/3A/005)) 

$$
{\bf W}_{\mathrm{p}} = 
\iota({\bf c}\otimes{\bf c}\cdot{\bf D} - {\bf D}\cdot{\bf c}\otimes{\bf c})
+
\zeta({\bf c}\otimes{\bf c}\cdot{\bf D}^2 - {\bf D}^2\cdot{\bf c}\otimes{\bf c})
.
$$

By requiring that basal planes preserve their orientation when subject to simple shear (like a deck of cards), it can be shown that $\iota=1$ and $\zeta=0$.

The corresponding effect on the continuous $c$-axis distribution is modelled as a conservative advection process on the surface of the unit sphere:

$$ \frac{\partial}{\partial t} n = -\nabla_{S^2}\cdot(n{\bf \dot{c}}) 
\qquad\Longrightarrow\qquad
\frac{\partial}{\partial t} {\boldsymbol  n} = {\bf M_{\mathrm{LROT}}} \cdot {\boldsymbol n} ,
$$

where ${\bf M_{\mathrm{LROT}}}$ is given analytically in [Rathmann et al. (2021)](https://doi.org/10.1017/jog.2020.117).

!!! tip "c-axis velocity field"
    The normalized $c$-axis velocity fields for the three modes of deformation considered are:
    
    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/latrot-caxis-velocity/latrot-caxis-velocity.png){: style="width:600px"}

### Example

```python
import numpy as np
from specfabpy import specfabpy as sf
# L=6 truncation is sufficient in this case, but larger L allows a very strong fabric to develop 
#  and minimizes the effect that regularization has on low wavenumber modes (l=2,4)
lm, nlm_len = sf.init(8) 

### Velocity gradient tensor experienced by parcel
ugrad = np.diag([0.5, 0.5, -1.0]) # Unconfined compression along z-direction (equal extension in x and y)
D = (ugrad+np.transpose(ugrad))/2 # Symmetric part (strain-rate tensor)
W = (ugrad-np.transpose(ugrad))/2 # Anti-symmetric part (spin tensor)

### Numerics 
Nt = 25 # Number of time steps
dt = 0.05 # Time-step size

### Initialize fabric as isotropic
nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # Array of expansion coefficients
nlm[0,0] = 1/np.sqrt(4*np.pi) # Normalized ODF at t=0

### Euler integration of lattice rotation + regularization
for tt in np.arange(1,Nt):
    nlm_prev = nlm[tt-1,:] # Previous solution
    iota, zeta = 1, 0 # "Deck of cards" behavior 
    M_LROT = sf.M_LROT(nlm_prev, D, W, iota, zeta) # Lattice rotation operator (nlm_len x nlm_len matrix) 
    M_REG  = sf.M_REG(nlm_prev, D)     # Regularization operator   (nlm_len x nlm_len matrix)
    M      = M_LROT + M_REG
    nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # Euler step
    nlm[tt,:] = sf.apply_bounds(nlm[tt,:]) # Apply spectral bounds if needed

# To plot the ODF at time-step tt, use plot_ODF(nlm[tt,:], lm, ax, geo) (see the corresponding demo)
# To derive the a2 structure tensor at time-step tt, use sf.a2(nlm[tt,:]) (see the corresponding demo)
```

See also [demo code in repository](https://github.com/nicholasmr/specfab/tree/main/demo/fabric-evolution).

- - -

## Discontinous dynamic recrystallization (DDRX) 

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/fabric-evolution/animation-DDRX.gif){: style="width:650px"}

DDRX is modelled as a grain orientation or mass decay/production process on the unit sphere ([Placidi and others, 2010](https://doi.org/10.1007/s00161-009-0126-0)):

$$ 
\frac{\partial}{\partial t} n = \Gamma n 
\qquad\Longrightarrow\qquad
\frac{\partial}{\partial t} {\boldsymbol  n} = {\bf M_{\mathrm{DDRX}}} \cdot {\boldsymbol n} ,
$$

where the decay rate 

$$\Gamma = \Gamma_0\left(D- {\langle} D {\rangle}\right) \quad\text{(decay rate)}$$

depends on the decay-rate magnitude $\Gamma_0$  (`Gamma0` in specfab), and the "deformability" $D$ as a function of the stress tensor ${\bf S}$:

$$
D = \frac{({\bf S}\cdot{\bf S}):({\bf c}\otimes{\bf c}) - {\bf S}:({\bf c}\otimes{\bf c}\otimes{\bf c}\otimes{\bf c}):{\bf S}}{{\bf S}:{\bf S}}\quad\text{(deformability)}
.
$$

!!! note 
    The average deformability, $\langle D\rangle$, depends on the instantaneous CPO state &mdash; specifically, the [structure tensors](cpo-representation.md) `a2` and `a4` &mdash; 
    making the corresponding matrix problem nonlinear by conserving the total number of grain orientations or mass density [depending on how normalization is interpreted](cpo-representation.md), the latter arguably resting on stronger physical grounds.

!!! tip "Decay/production rate"

    The normalized DDRX decay rate $\Gamma/\Gamma_0 = D - \langle D \rangle$ is an orientation dependent field that favors nucleation (orientation/mass production) in the directions where the resolved basal-plane shear stress is maximal, and orientation/mass decay elsewhere. 

    The the normalized decay rate for the three modes of deformation considered are:

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/ddrx-decayrate/ddrx-decayrate.png){: style="width:600px"}

### Example

```python
import numpy as np
from specfabpy import specfabpy as sf
# L=6 truncation is sufficient in this case, but larger L allows a very strong fabric to develop 
#  and minimizes the effect that regularization has on low wavenumber modes (l=2,4)
lm, nlm_len = sf.init(8) 

### Stress tensor experienced by parcel
S = np.diag([0.0, 1.0, -1.0]) # Confined compression along z-direction (extension confined to y-direction)
#S = np.diag([0.5, 0.5, -1.0]) # Unconfined compression along z-direction
Gamma0 = 10 # DDRX decay-rate magnitude (may depend on temperature, strain-rate, and other factors, see e.g. Richards et al., 2021)

### Numerics 
Nt = 25 # Number of time steps
dt = 0.05 # Time-step size

### Initialize fabric as isotropic
nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # Array of expansion coefficients
nlm[0,0] = 1/np.sqrt(4*np.pi) # Normalized ODF at t=0

### Euler integration of DDRX
for tt in np.arange(1,Nt):
    nlm_prev = nlm[tt-1,:] # Previous solution
    M = Gamma0 * sf.M_DDRX(nlm_prev, S) # DDRX operator (nlm_len x nlm_len matrix)
    nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # Complete Euler step
    nlm[tt,:] = sf.apply_bounds(nlm[tt,:]) # Apply spectral bounds if needed    

# To plot the ODF at time-step tt, use plot_ODF(nlm[tt,:], lm, ax, geo) (see the corresponding demo)
# To derive the a2 structure tensor at time-step tt, use sf.a2(nlm[tt,:]) (see the corresponding demo)
```

- - -

## Continous dynamic recrystallization (CDRX) 

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/fabric-evolution/animation-CDRX.gif){: style="width:650px"}

Polygonization (rotation recrystallization, CDRX) accounts for the division of grains along internal sub-grain boundaries when exposed to bending stresses. 
In effect, CDRX reduces the average grain size upon grain division but does not necessarily change the CPO much ([Alley, 1992](https://doi.org/10.3189/S0022143000003658)). 

The model follows [Gödert (2003)](https://doi.org/10.1007/s001610050095) by approximating the effect of CDRX on the distribution of grain orientations/mass as a Laplacian diffusive process on $S^2$:

$$
\frac{\partial}{\partial t} n = \Lambda\nabla^2 n 
\qquad\Longrightarrow\qquad
\frac{\partial}{\partial t} {\boldsymbol  n} = {\bf M_{\mathrm{CDRX}}} \cdot {\boldsymbol n} .
$$

### Example 

To model CDRX, add the following contribution to the total fabric operator ${\bf M}$:
```python
M += Lambda*sf.M_CDRX(nlm)
```
where `Lambda` ($\Lambda$) is the CDRX rate-factor magnitude that possibly depends on temperature, stress, strain-rate, etc. ([Richards et al., 2021](https://doi.org/10.1016/j.epsl.2020.116718)).

- - -

## Regularization

![](https://raw.githubusercontent.com/nicholasmr/specfab/8afe59b8847e16761d69abcc5c8ec0327c85e61c/tests/calibrate-regularization/animation.gif){: style="width:700px"}

As $n(\theta,\phi)$ becomes anisotropic due to CPO processes, the coefficients $n_l^m$ associated with high wavenumber modes (large $l$ and $m$, and thus small-scale structure) must increase in magnitude relative to the low wavenumber coefficients (small $l$ and $m$). 

One way to visualize this is by the [angular power spectrum](https://en.wikipedia.org/wiki/Spherical_harmonics#Spectrum_analysis)

$$ 
S(l) = \frac{1}{2l + 1} \sum_{m=-l}^l \left\vert n_l^m \right\vert^2 ,
$$

which grows with time. 
In the animation above, the left-hand panel shows how the power spectrum evolves under lattice rotation (unconfined vertical compression) compared to the end-member case of a delta-function (dashed line).

If the expansion series is truncated at $l=L$, then $l{\gt}L$ modes cannot evolve, and the truncated solution will reach an unphysical quasi-steady state. To prevent this, regularization must be introduced.
Specfab uses Laplacian hyper diffusion ($k>1$) as regularization in $S^2$

$$ 
\frac{\partial}{\partial t} n_l^m={\nu}[l(l+1)]^{k} n_l^m 
\qquad\Longrightarrow\qquad
\frac{\partial}{\partial t} {\boldsymbol  n} = {\bf M_{\mathrm{REG}}} \cdot {\boldsymbol n} 
,
$$

that can be added to the fabric evolution operator ${\bf M}$ as follows:
```python
M += sf.M_REG(nlm, D)
```
This allows the growth of high wavenumber modes to be disproportionately damped (green line compared to red line in animation above).
However, `sf.M_REG()` is currently only calibrated for $L = 4,6,8,20$ &mdash; that is, optimal values of $\nu$ and $k$ are only available for these truncations.

!!! note 
    As a rule-of-thumb, regularization affects the highest and next-highest modes $l{\geq}L-2$ and can therefore not be expected to evolve freely. 
    This, in turn, means that [structure tensors](cpo-representation.md) `a2` and `a4`, and hence [calculated enhancement factors](enhancements-strainrate.md), will be affected by regularization unless $L{\geq}8$ is chosen. 


