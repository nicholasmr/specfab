# Strain-rate enhancements

Given an anisotropic rheology ${\bf D}({\bf S})$, where ${\bf D}$ and ${\bf S}$ are the 
strain-rate and deviatoric stress tensors, respectively, 
the *directional strain-rate enhancement factors* $E_{ij}$ are defined as the $({\bf e}_i, {\bf e}_j)$-components of ${\bf D}$ relative to that of the rheology in the isotropic limit (isotropic CPO):

$$ 
E_{ij} = \frac{
{\bf e}_i \cdot {\bf D}({\bf S}) \cdot {\bf e}_j 
}{
{\bf e}_i \cdot {\bf D}_{\mathrm{iso}}({\bf S}) \cdot {\bf e}_j 
}
,
 \qquad(1)
$$

for a stress state aligned with $({\bf e}_i, {\bf e}_j)$:

$$
{\bf S}({\bf e}_i, {\bf e}_j) = \tau_0
\begin{cases}
    {\bf I} - 3{\bf e}_i \otimes {\bf e}_i \;\;\quad\quad\text{if}\quad i=j \\
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

## Eigenenhancements

*Eigenenhancements* are defined as the enhancement factors w.r.t. the CPO symmetry axes (${\bf m}_i$): 

$${\bf e}_i = {\bf m}_i .$$

These are the enhancements factors needed to specify the viscous anisotropy in [bulk rheologies](constitutive-viscoplastic.md):

| Transversely isotropic | Orthotropic |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous.png){: style="width:260px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous.png){: style="width:350px"} |

## Grain homogenization 

Calculating $E_{ij}$ using (1) for a given CPO requires an *effective* rheology that takes the microstructure into account.

In the simplest case, polycrystals may be regarded as an ensemble of interactionless grains (monocrystals), subject to either a homogeneous stress field over the polycrystal scale:

$$
{\bf S}' = {\bf S}
,
\qquad\qquad \text{(Sachs's hypothesis)}
$$

or a homogeneous stain-rate field:

$$
{\bf D}' = {\bf D} 
,
\qquad\qquad \text{(Taylor's hypothesis)}
$$

where ${\bf S}'$ and ${\bf D}'$ are the *microscopic* (grain-scale) stress and strain-rate tensors, respectively.

The effective rheology can then be approximated as the ensemble-averaged monocrystal rheology for either case:

$$
{\bf D}^{\mathrm{Sachs}} = \langle {\bf D}'({\bf S}') \rangle = \langle {\bf D}'({\bf S}) \rangle
,
\qquad\qquad \text{(Sachs homogenization)}
$$

$$
\qquad
{\bf D}^{\mathrm{Taylor}} = \langle {\bf S}'({\bf D}') \rangle^{-1} = \langle {\bf S}'({\bf D}) \rangle^{-1}
,
\qquad \text{(Taylor homogenization)}
$$

where $\langle \cdot \rangle^{-1}$ inverts the tensorial relationship.

If a linear combination of the two homogenizations is considered, equation (1) can be written as 

$$
E_{ij} = (1-\alpha) \frac{{\bf e}_i \cdot {\bf D}^{\mathrm{Sachs}}({\bf S}) \cdot {\bf e}_j}
{{\bf e}_i \cdot {\bf D}^{\mathrm{Sachs}}_{\mathrm{iso}}({\bf S}) \cdot {\bf e}_j }
+ {\alpha} \frac{ {\bf e}_i \cdot {\bf D}^{\mathrm{Taylor}}({\bf S}) \cdot {\bf e}_j }
{ {\bf e}_i \cdot {\bf D}^{\mathrm{Taylor}}_{\mathrm{iso}}({\bf S}) \cdot {\bf e}_j }
,
$$

or simply

$$
E_{ij} = (1-\alpha)E_{ij}^{\mathrm{Sachs}} + {\alpha}E_{ij}^{\mathrm{Taylor}} ,
$$

where $\alpha$ is a free parameter.

!!! warning "Grain parameters"
    The grain viscous parameters used for homogenization should be understood as the *effective* values needed to reproduce deformation experiments on polycrystals; they are not the values derived from experiments on single crystals.

### Transversely isotropic grains

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous-monocrystal.png){: style="width:210px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polycrystal.png){: style="width:220px"} |

If grains are approximately transversely isotropic, the grain rheology can be modelled using the [transversely isotropic power-law rheology](constitutive-viscoplastic.md).
This requires specifying the grain eigenenhancements $E_{mm}'$ and $E_{mt}'$, the power-law exponent $n'$, and the Taylor&mdash;Sachs weight $\alpha$.

#### Example for glacier ice

For glacier ice, we follow the literature and rename 

$$
{\bf c} = {\bf m}^\prime, \\
{\bf a} = {\bf t}^\prime.
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

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/enhancement-factor/ice/Eij-state-space/Eij-state-space.gif){: style="width:750px"}

    Clearly, there is a tradeoff between how shear enhanced ($E_{mt}$) and how hard for axial compression ($E_{mm}$) the model allows a unidirectional CPO to be.
    

!!! tip "Evolving CPO"

    The below animation shows the directional enhancement factors for a CPO evolving under [uniaxial compression](deformation-modes.md) along ${\hat {\bf z}}$ when subject to [lattice rotation](cpo-dynamics-tranisotropic.md).
    Enhancement factors are calculated w.r.t. the spherical coordinate basis vectors $({\bf e}_1, {\bf e}_2, {\bf e}_3) = ({\hat{\bf r}},{\hat{\boldsymbol \theta}},{\hat{\boldsymbol \phi}})$.

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/cube-crush-animation/S2-maps/S2-Eij.gif){: style="width:660px"}

### Orthotropic grains

| Monocrystal | Polycrystal |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous-monocrystal.png){: style="width:250px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/polycrystal.png){: style="width:220px"} |

If grains are approximately orthotropic, the grain rheology can be modelled using the [orthotropic power-law rheology](constitutive-viscoplastic.md).
This requires specifying the grain eigenenhancements $E_{ij}'$, the power-law exponent $n'$, and the Taylor&mdash;Sachs weight $\alpha$.

#### Example for olivine

Not yet available.

