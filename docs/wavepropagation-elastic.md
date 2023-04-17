# Elastic wave propagation

Elastic P and S plan-wave velocities, permitted in a polycrystal with an arbitrary CPO, can be determined for any linear combination of the Voigt and Reuss homogenization schemes by solving for the eigenvalues of the [acoustic tensor](https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Elasticity/Elasticity.htm):

$$
{\bf {Q}}= (1-\alpha){\bf Q_{\mathrm{Reuss}}} + \alpha{\bf Q_{\mathrm{Voigt}}},
$$
which depends on the CPO, grain elastic parameters, and direction of propagation.

!!! note "Voigt and Reuss homogenizations"
    The Voigt scheme ($\alpha=1$) assumes the strain field is homogeneous over the polycrystal scale, whereas the Reuss scheme ($\alpha=0$) assumes the stress is homogeneous.
    Moreover, in these homogenizations, grains are assumed interactionless and the bulk elastic behaviour is therefore simply the grain-orientation-averaged elastic behaviour subject to either homogeneous stress or strain assumptions over the polycrystal scale.

!!! warning "Grain parameters" 
    The grain elastic parameters should be understood as the *effective* polycrystal values needed to reproduce experimental results, and not measured values derived from experiments on single crystals.

## Transversely isotropic grains

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-elastic-monocrystal.png){: style="width:210px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polycrystal.png){: style="width:220px"} |

If grains are approximately transversely isotropic, the grain elastic behaviour can be modelled using the [transversely isotropic elastic constitutive equation](constitutive-elastic.md).
This requires specifying the three grain parameters $\lambda'$, $\mu'$, $\hat{\lambda}'$, $\hat{\mu}'$, $\hat{\gamma}'$, and the Voigt&mdash;Reuss weight $\alpha$.

### Example for glacier ice

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(4) # L=4 truncation is sufficient in this case
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients

# c-axis number distribution (nlm) from fourth-order structure tensor (a4)
p = np.array([0,0,1]) # preferred c-axis direction
a4 = np.einsum('i,j,k,l', p,p,p,p) # a4 if ODF = deltafunc(r-p) 
nlm[:] = sf.a4_to_nlm(a4) # l=2,4 expansion coefficients for corresponding ODF (a4 is normalized)

# Physical parameters
rho = 917 # density of ice
C11,C33,C55,C12,C13 = 14.060e9, 15.240e9, 3.060e9, 7.150e9, 5.880e9 # Bennett (1968) parameters
lam,mu,Elam,Emu,Egam = sf.Cij_to_Lame_tranisotropic(C11,C33,C55,C12,C13) 

# Homogenization scheme
alpha = 0.5 # Voigt--Reuss weight, where 0.5 = Hill average

# Propagation directions of interest
theta, phi = np.deg2rad([90,70,]), np.deg2rad([0,10,]) # wave-vector directions (theta is colatitude, phi is longitude)

# Calculate phase velocities
Vi = sf.Vi_elastic_tranisotropic(nlm, alpha, lam,mu,Elam,Emu,Egam, rho, theta,phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
```

## Orthotropic grains

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-elastic-monocrystal.png){: style="width:250px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/polycrystal.png){: style="width:220px"} |

If grains are approximately orthotropic, the grain elastic behaviour can be modelled using the [orthotropic elastic constitutive equation](constitutive-elastic.md).
This requires specifying the ten grain parameters $\lambda_{11}'$, $\lambda_{22}'$, $\lambda_{33}'$, $\lambda_{12}'$, $\lambda_{13}'$, $\lambda_{23}'$, $\mu_{1}'$, $\mu_{2}'$, $\mu_{3}'$, and the Voigt&mdash;Reuss weight $\alpha$.

### Example for olivine

Not yet available.

<!--
For polycrystals with orthotropic grains characterized by *effective* values of the grain elastic constants $a,b,c$ (see [elasticities](constitutive-elastic.md)).

**Example for olivine**

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(4) # L=4 truncation is sufficient in this case
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients

# c-axis number distribution (nlm) from fourth-order structure tensor (a4)
p = np.array([0,0,1]) # preferred c-axis direction
a4 = np.einsum('i,j,k,l', p,p,p,p) # a4 if ODF = deltafunc(r-p) 
nlm[:] = sf.a4_to_nlm(a4) # l=2,4 expansion coefficients for corresponding ODF (a4 is normalized)

# Physical parameters
rho = 917 # density of ice
C11,C33,C55,C12,C13 = 14.060e9, 15.240e9, 3.060e9, 7.150e9, 5.880e9 # Bennett (1968) parameters
lam,mu,Elam,Emu,Egam = sf.Cij_to_Lame_tranisotropic(C11,C33,C55,C12,C13) 

# Homogenization scheme
alpha = 0.5 # Voigt--Reuss weight, where 0.5 = Hill average

# Phase velocities
theta, phi = np.deg2rad([90,70,]), np.deg2rad([0,10,]) # wave-vector directions (theta is colatitude, phi is longitude)
Vi = sf.Vi_elastic_tranisotropic(nlm, alpha, lam,mu,Elam,Emu,Egam, rho, theta,phi) # phase velocities are V_S1=vi[0,:], V_S2=vi[1,:], V_P=vi[2,:]
```
-->
