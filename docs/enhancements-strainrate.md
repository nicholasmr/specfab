# Strain-rate enhancements

## Definition

Given an anisotropic rheology ${\bf D}({\bf S})$, where ${\bf D}$ and ${\bf S}$ are the 
strain-rate and deviatoric stress tensors, respectively, 
the *directional strain-rate enhancement factors*, $E_{ij}$, are defined as the $({\bf e}_i, {\bf e}_j)$-components of ${\bf D}$ relative to that of the rheology in the isotropic limit (isotropic CPO):

$$ 
E_{ij} = \frac{
{\bf e}_i \cdot {\bf D}({\bf S}({\bf e}_i, {\bf e}_j)) \cdot {\bf e}_j 
}{
{\bf e}_i \cdot {\bf D}_{\mathrm{iso}}({\bf S}({\bf e}_i, {\bf e}_j)) \cdot {\bf e}_j 
}
,
 \qquad(1)
$$

for a stress state aligned with $({\bf e}_i, {\bf e}_j)$:

$$
{\bf S}({\bf e}_i, {\bf e}_j) = \tau_0
\begin{cases}
    {\bf I}/3 - {\bf e}_i \otimes {\bf e}_i \;\;\quad\quad\text{if}\quad i=j \\
    {\bf e}_i \otimes {\bf e}_j + {\bf e}_j \otimes {\bf e}_i \quad\text{if}\quad i\neq j \\
\end{cases}
.
$$

In this way:

* ${E_{11}}$ is the longitudinal strain-rate enhancement along ${\bf e}_{1}$ when subject to compression along ${\bf e}_{1}$

* ${E_{12}}$ is the ${\bf e}_{1}$&mdash;${\bf e}_{2}$ shear strain-rate enhancement when subject to shear in the ${\bf e}_{1}$&mdash;${\bf e}_{2}$ plane

and so on.

!!! note "Hard or soft"

    $E_{ij}>1$ implies the material response is *softened* due to fabric (compared to an isotropic CPO), whereas $E_{ij}<1$ implies *hardening*.

### Eigenenhancements

*Eigenenhancements* are defined as the enhancement factors w.r.t. the CPO symmetry axes (${\bf m}_i$): 

$${\bf e}_i = {\bf m}_i .$$

These are the enhancements factors needed to specify the viscous anisotropy in [bulk rheologies](constitutive-viscoplastic.md):

| Transversely isotropic | Orthotropic |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous.png){: style="width:260px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous.png){: style="width:350px"} |

## Grain homogenization schemes 

Using (1) to calculate $E_{ij}$ for a given CPO requires an effective rheology that takes the microstructure into account.

In the simplest case, polycrystals may be regarded as an ensemble of interactionless grains (monocrystals) subject to either a homogeneous strain field (Taylor's hypothesis) or homogeneous stress field (Sachs's hypothesis) over the polycrystal scale. 
In this way, the effective rheology is simply the grain-orentation-averaged rheology, assuming homogeneous stress or strain-rate over the polycrystal scale.

Any linear combination of the two homogenizations is supported:

$$ 
E_{ij} = (1-\alpha)E_{ij}^{\mathrm{Sachs}} + {\alpha}E_{ij}^{\mathrm{Taylor}} ,
$$

where $E_{ij}^{\mathrm{Sachs}}$ and $E_{ij}^{\mathrm{Taylor}}$ are calculated with (1) assuming constant $\bf{S}$ and $\bf{D}$, respectively.

!!! warning "Grain parameters"
    The grain viscous parameters below should be understood as the *effective* polycrystal values needed to reproduce deformation experiments, and not measured values derived from experiments on single crystals.

### Transversely isotropic grains

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous-monocrystal.png){: style="width:210px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polycrystal.png){: style="width:220px"} |

If grains are approximately transversely isotropic, the grain rheology can be modelled using the [transversely isotropic power-law rheology](constitutive-viscoplastic.md).
This requires specifying the three grain parameters $n$, $E_{cc}$ and $E_{ca}$, and the Taylor&mdash;Sachs weight $\alpha$.

#### Example for glacier ice

The below example for glacier ice shows how $E_{ij}$ may be calculated given `a2`, `a4`, or `nlm`. 

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(4) # L=4 truncation is sufficient in this case

### Make synthetic ODF
# Unidirectional CPO: all c-axes aligned in z-direction (ODF = deltafunc(r-m))
m = np.array([0,0,1]) 
nlm = np.zeros((nlm_len), dtype=np.complex64) # Array of expansion coefficients
if True: # use a2
    a2 = np.einsum('i,j', m,m) # Outer product
    nlm[:6] = sf.a2_to_nlm(a2) # Derive corresponding expansion coefficients
else: # use a4
    a4 = np.einsum('i,j,k,l', m,m,m,m) # Outer product
    nlm[:15] = sf.a4_to_nlm(a4) # Derive corresponding expansion coefficients

### Coordinate basis vectors for enhancement-factor calculations
(e1,e2,e3, eigvals) = sf.frame(nlm, 'e') # a2 eigen basis
#(e1,e2,e3) = (np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])) # x,y,z cartesian basis

### Transversely isotropic monocrystal parameters for ice (Rathmann & Lilien, 2021)
n     = 1      # n=1 => linear grain rheology (nonlinear not fully supported)
Ecc   = 1      # Enhancement factor for compression along c-axis
Eca   = 1e3    # Enhancement factor for shear parallel to basal plane
alpha = 0.0125 # Taylor--Sachs weight

### Calculate enhancement-factor matrix in the basis (e1,e2,e3)
Eij = sf.Eeiej(nlm, e1,e2,e3, Ecc,Eca,alpha,n) 
```

### Orthotropic grains

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous-monocrystal.png){: style="width:250px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/polycrystal.png){: style="width:220px"} |

If grains are approximately orthotropic, the grain rheology can be modelled using the [orthotropic power-law rheology](constitutive-viscoplastic.md).
This requires specifying the eight grain parameters $n$, $E_{11}$, $E_{22}$, $E_{33}$, $E_{12}$, $E_{13}$, $E_{23}$, and the Taylor&mdash;Sachs weight $\alpha$.

#### Example for olivine

Not yet available.

