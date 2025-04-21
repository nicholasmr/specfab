# Viscous anisotropy of glacier ice 

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous-monocrystal.png){: style="width:120px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polycrystal.png){: style="width:200px"} |

If grains are approximately transversely isotropic, the grain rheology can be modelled using the [transversely isotropic power-law rheology](constitutive-viscoplastic.md).
This requires specifying the grain eigenenhancements $E_{mm}'$ and $E_{mt}'$, the power-law exponent $n'$, and the Taylor&mdash;Sachs weight $\alpha$.

For glacier ice, we follow the literature and rename 

$$
{\bf c} = {\bf m}^\prime \quad\text{and}\quad {\bf a} = {\bf t}^\prime.
$$

The below code example shows how to calculate $E_{ij}$ given `a4` (or `nlm`) assuming the grain parameters proposed by [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88): 

```python
import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(8) 

### Synthetic unidirectional CPO (all c-axes aligned in z-direction)
m = np.array([0,0,1]) 
a4 = np.einsum('i,j,k,l', m,m,m,m) # 4-times repeated outer product of m
nlm = np.zeros((nlm_len), dtype=np.complex64)
nlm[:sf.L4len] = sf.a4_to_nlm(a4) # derive corresponding expansion coefficients

### Basis for enhancement factor calculations
(e1,e2,e3, eigvals) = sf.frame(nlm, 'e') # enhancement factors are w.r.t. a^(2) basis (i.e. eigenenhancements)
#(e1,e2,e3) = np.eye(3) # enhancement factors are w.r.t. Cartesian basis (x,y,z)

### Transversely isotropic monocrystal parameters for ice (Rathmann & Lilien, 2021)
n_grain   = 1        # power-law exponent: n=1 => linear grain rheology, nonlinear (n>1) is unsupported.
Eij_grain = (1, 1e3) # grain eigenenhancements (Ecc,Eca) for compression along c-axis (Ecc) and for shear parallel to basal plane (Eca)
alpha     = 0.0125   # Taylor--Sachs weight

### Calculate enhancement factors w.r.t. (e1,e2,e3)
Eij = sf.Eij_tranisotropic(nlm, e1,e2,e3, Eij_grain,alpha,n_grain) # Eij=(E11,E22,E33,E23,E13,E12)
```

!!! tip "Choosing grain parameters for glacier ice"

    The grain parameters proposed by [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88) assume a linear viscous ($n'=1$) response and promote the activation of basal glide by making that slip system soft compared to other systems: $E_{ca}' > 1$, whereas $E_{cc}'=1$. 
    This reduces the problem to that of picking $E_{ca}'$ and $\alpha$, which [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88) chose such that deformation tests on unidirectional CPOs (perfect single maximum) are approximately reproduced: $E_{mt}=10$ while $E_{mt}/E_{pq} \sim 10^4$, where $p,q$ denote directions at $45^\circ$ to ${\bf m}$.

    The effect of choosing alternative $E_{ca}'$ and $\alpha$ (left plot) on eigenenhancements for different CPO states (right plot) is here shown for combinations of $E_{ca}'$ and $\alpha$ that fulfill $E_{mt}=10$ for a unidirectional CPO:

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/research/calibrate-Eij/ice/Eij-state-space/Eij-state-space.gif){: style="width:750px"}

    Clearly, there is a tradeoff between how shear enhanced ($E_{mt}$) and how hard for axial compression ($E_{mm}$) the model allows a unidirectional CPO to be.
    


