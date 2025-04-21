## Discontinous dynamic recrystallization (DDRX) 

Following [Placidi and others (2010)](https://doi.org/10.1007/s00161-009-0126-0), DDRX is modeled as a spontaneous mass decay&mdash;production process in orientation space $S^2$, intended to represent the combined effect of nucleation and grain boundary migration. 
That is, mass is spontaneously exchanged between grains with different orientations depending on the local stress state, strain rate, and temperature, in a statistical sense. 

The decay&mdash;production rate is defined as

$$
\Gamma = \Gamma_0\left(D- {\langle} D {\rangle}\right) 
$$

where the  prefactor $\Gamma_0$ accounts for the preferential (Arrhenius) activation at warm temperatures and the effect of strain-rate magnitude, defined as
\begin{align}
\Gamma_0 = \sqrt{\frac{{\bf D}:{\bf D}}{2}} A_{\Gamma}\exp(-Q_{\Gamma}/RT)
.
\end{align}
Here, ${\bf D}$ is the strain-rate tensor, $R$ is the gas constant, $T$ is the temperature, and $A_{\Gamma}$ and $Q_{\Gamma}$ have been calibration by [Richards et al. (2021)](https://doi.org/10.1016/j.epsl.2020.116718) and [Lilien et al. (2023)](https://doi.org/10.1017/jog.2023.78). 

The deformability $D$ (not to be mistaken for the strain-rate tensor ${\bf D}$) is the normalized square of the basal-plane resolved shear stress

$$
D = 5\frac{({\bf S}\cdot{\bf S}):({\bf n}\otimes{\bf n}) - {\bf S}:({\bf n}\otimes{\bf n}\otimes{\bf n}\otimes{\bf n}):{\bf S}}{{\bf S}:{\bf S}},
$$

where ${\bf n}$ is an arbitrary slip-plane normal ($c$-axis for ice) and the factor of 5 is conventional to include.
Because $D$ is largest for grains with an orientation favorable to easy glide, mass is spontaneously created/added to grains with such preferred orientations (in a statistical sense).
Conversely, mass spontaneously decays if $D<{\langle} D {\rangle}$, corresponding to grains with an unfavorable orientation being consumed by grains with a more favorable orientation to basal glide. 
Here, ${\langle} D {\rangle}$ is the grain-average deformability of the polycrystal. 

!!! warning "Nonlinear process"

    Since average deformability $\langle D\rangle$ depends on the instantaneous CPO state &mdash; specifically, the [structure tensors](cpo-representation.md) `a2` and `a4` &mdash; this crystal process is nonlinear (renders a nonlinear matrix problem below).

!!! warning "Limited use"

    Notice that this model is currently only relevant to slip-system normals. 
    The model is therefore not yet useful for e.g. olivine. 

!!! note "Glacier ice"

    The normalized decay&mdash;production rate is shown below for three different stress states:

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation-modes/stress-modes.png#center){: style="width:550px"}

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/ddrx-decayrate/ddrx-decayrate.png){: style="width:560px"}


### Matrix model 

The corresponding effect on the continuous distribution function is 

$$ 
\frac{\mathrm{D} n}{\mathrm{D} t} = \Gamma n 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}_n}{\mathrm{D} t} = {\bf M_{\mathrm{DDRX}}} \cdot {\bf s}_n ,
$$

where ${\bf M_{\mathrm{DDRX}}}$ is given analytically in [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88).

### Code example

```python
import numpy as np
from specfabpy import specfab as sf
# L=8 truncation is sufficient in this case, but larger L allows a very strong fabric to  
#  develop and minimizes the effect that regularization has on low wavenumber modes (l=2,4)
lm, nlm_len = sf.init(8) 

### Stress tensor experienced by parcel

S = np.diag([0.5, 0.5, -1.0]) # uniaxial compression along z-axis

### Numerics 

Nt = 25   # number of time steps
dt = 0.05 # time-step size

### Initial fabric state

nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # n state vector
nlm[0,0] = 1/np.sqrt(4*np.pi) # normalized isotropic state at t=0

### Euler integration

for tt in np.arange(1,Nt):
    nlm_prev = nlm[tt-1,:] # previous solution
    Gamma0 = 10            # DDRX decay-rate magnitude
    M = Gamma0 * sf.M_DDRX(nlm_prev, S) # DDRX operator
    nlm[tt,:] = nlm_prev + dt*np.matmul(M, nlm_prev) # euler step
    nlm[tt,:] = sf.apply_bounds(nlm[tt,:]) # apply spectral bounds if needed    
```

