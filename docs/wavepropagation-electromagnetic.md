# Electromagnetic wave propagation

## Problem

We seek plane wave solutions to Maxwell's equations in a non-conducting, source-free, anisotropic linear dielectric medium

$$
\nabla^2 {\bf E} = \mu {\boldsymbol \epsilon} \frac{\partial^2 {\bf E}}{\partial t^2},
$$

where ${\bf E}$ is the electric field, ${\boldsymbol \epsilon}$ is the bulk dielectric permittivity tensor, and $\mu$ the bulk isotropic permeability of the medium.

Substituting ${\bf E}$ for a plane wave solution, ${\bf E} = {\bf E}_0 \exp[i({\bf k}\cdot {\bf x} - \omega t)]$, the problem reduces to

$$
({\bf K} + \omega^2 \mu {\boldsymbol \epsilon}) {\bf E} = {\bf 0},
$$

where ${\bf K} = {\bf k}\times{\bf k} \times$ is the [matrix representation](https://en.wikipedia.org/wiki/Cross_product#Alternative_ways_to_compute) of the twice-applied cross product.

The above equation requires 

$$
\det( {\bf K} + \omega^2 \mu {\boldsymbol \epsilon} ) = 0
.
$$
<!--{\boldsymbol \epsilon}^{-1} \frac{\hat{{\bf k}}^2}{\mu} {\bf E} = \frac{\omega^2}{k^2} {\bf E} -->

Evidently, the eigenvalues and eigenvectors of ${\boldsymbol \epsilon}^{-1} {\bf K}/(\mu k^2)$ are the permitted wave velocities squared and wave polarization, respectively, where

$$
V^2 = \frac{\omega^2}{k^2}
$$

is the wave velocity squared.

!!! note "Homogenization"

    If wave lengths are much longer than the average grain size, the problem may be close by approximating the bulk permittivity tensor of a polycrystal by the grain-average permittivity tensor

    $${\boldsymbol \epsilon} = \langle {\boldsymbol \epsilon}' \rangle,$$
    
    constructed by averaging over all grain orientations (over the CPO).

## Transversely isotropic grains

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-electromagnetic-monocrystal.png){: style="width:210px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polycrystal.png){: style="width:220px"} |


If grains are approximately transversely isotropic w.r.t. the symmetry axis ${\bf m}'$, the dielectric permittivity tensor of a single crystal can be written as
$$
{\boldsymbol \epsilon}' =  (2\epsilon_{t}' + \epsilon_{m}') \frac{\bf I}{3}
+ (\epsilon_{m}'-\epsilon_{t}') \left( {\bf m}'^2 - \frac{\bf I}{3} \right),
$$
where $\epsilon_{m}'$ and $\epsilon_{t}'$ are the permittivities parallel and perpendicular to the symmetry axis.
In this case, the grain-averaged permitivity is simply 

$$
{\boldsymbol \epsilon} = (2\epsilon_{t}' + \epsilon_{m}') \frac{\bf I}{3}
+ (\epsilon_{m}'-\epsilon_{t}') \left( \langle {\bf m}'^2 \rangle - \frac{\bf I}{3} \right)
,
$$

where $\langle {\bf m}'^2 \rangle$ is the [second-order structure tensor](cpo-structuretensors.md) (aka ${\bf a}^{(2)}$).

### Example for glacier ice

To be documented.

<!--

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(2) # L=2 is sufficient here
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients

# Physical parameters
# ... note eps=eps0*epsr, mu=mu0*mur
epsr_grain = (3.17, 3.17-0.034) # (epsc, epsa); relative permittivity of a single grain parallel (c) and perpendicular (a) to symmetry axis
mur = 1 # relative permeability of a single grain

# c-axis number distribution (nlm) from second-order structure tensor (a2)
p = np.array([0,0,1]) # preferred c-axis direction
a2 = np.einsum('i,j', p,p) # a2 if ODF = deltafunc(r-p) 
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients for corresponding normalized ODF

# Propagation directions of interest
theta, phi = np.deg2rad([90,70,]), np.deg2rad([0,10,]) # wave-vector directions (theta is colatitude, phi is longitude)

# Calculate phase velocities
Vi = sf.Vi_electromagnetic_tranisotropic(nlm, eps_grain, mu, theta,phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:]
```
-->

## Orthotropic grains

Not yet supported.
