# Strain rate enhancement

Given a bulk anisotropic rheology ${\bf D}({\bf S})$, where ${\bf D}$ and ${\bf S}$ are the 
bulk strain-rate and deviatoric stress tensors, respectively, 
the *directional strain-rate enhancement factors* $E_{ij}$ are defined as the $({\bf e}_i, {\bf e}_j)$-components of ${\bf D}$ relative to that of the rheology in the limit of an isotropic CPO:

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

## Eigenenhancements

*Eigenenhancements* are defined as the enhancement factors w.r.t. the CPO symmetry axes (${\bf m}_i$): 

$${\bf e}_i = {\bf m}_i .$$

These are the enhancements factors needed to specify the viscous anisotropy in [bulk rheologies](constitutive-viscoplastic.md):

| Transversely isotropic | Orthotropic |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous.png){: style="width:260px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous.png){: style="width:350px"} |

## Grain homogenization 

### Taylor&mdash;Sachs

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

If a linear combination of the two homogenizations is considered, equation (1) can be approximated as 

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

### Azuma&mdash;Placidi

🚧 *Not yet documented.*
