# Strain rate enhancement

Given a bulk anisotropic rheology ${\dot{\boldsymbol\epsilon}}(\boldsymbol\tau)$, where ${\dot{\boldsymbol\epsilon}}$ and ${\boldsymbol\tau}$ are the 
bulk strain-rate and deviatoric stress tensors, respectively, 
the *directional strain-rate enhancement factors* $E_{ij}$ are defined as the $({\bf e}_i, {\bf e}_j)$-components of ${\dot{\boldsymbol\epsilon}}$ relative to that of the rheology in the limit of an isotropic CPO:

$$ 
E_{ij} = \frac{
{\bf e}_i \cdot {\dot{\boldsymbol\epsilon}}({\boldsymbol\tau}) \cdot {\bf e}_j 
}{
{\bf e}_i \cdot {\dot{\boldsymbol\epsilon}}_{\mathrm{iso}}({\boldsymbol\tau}) \cdot {\bf e}_j 
}
,
 \qquad(1)
$$

for a stress state aligned with $({\bf e}_i, {\bf e}_j)$:

$$
{\boldsymbol\tau}({\bf e}_i, {\bf e}_j) = \tau_0
\begin{cases}
    {\bf I}/3 - {\bf e}_i {\bf e}_i \;\quad\text{if}\quad i=j \\
    {\bf e}_i {\bf e}_j + {\bf e}_j {\bf e}_i \quad\text{if}\quad i\neq j \\
\end{cases}
.
$$

In this way:

* ${E_{11}}$ is the longitudinal strain-rate enhancement along ${\bf e}_{1}$ when subject to compression along ${\bf e}_{1}$,

* ${E_{12}}$ is the ${\bf e}_{1}$&mdash;${\bf e}_{2}$ shear strain-rate enhancement when subject to shear in the ${\bf e}_{1}$&mdash;${\bf e}_{2}$ plane,

* ...and so on.

To be clear, $E_{ij}>1$ implies the material response is *softened* due to fabric (compared to an isotropic CPO), whereas $E_{ij}<1$ implies *hardening*.

!!! note "Example: Glacier ice"

    In the case of glacier ice with a strongly-developed preferred $c$-axis direction (*single maximum*; left figure below), $E_{ij}$ have been measured in lab tests of compression and shear along the preferred direction (right figure below): 

    /// html | div[style='float: left; width: 50%;']
    <br>
    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/deck-of-cards/ice-cards-2.png){: style="width:370px"}
    ///

    /// html | div[style='float: right; width: 44%;']
    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/deck-of-cards/ice-cards-VA.png){: style="width:350px"}
    ///

    /// html | div[style='clear: both;']
    ///

## Eigenenhancements

*Eigenenhancements* are defined as the enhancement factors with respect to the symmetry axes ${\bf m}_i$ of the CPO: 

$${\bf e}_i = {\bf m}_i .$$

These are the enhancements factors needed to specify the viscous anisotropy in [bulk rheologies](constitutive-viscoplastic.md):

| Transversely isotropic | Orthotropic |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/material-symmetries/icesym-traniso-viscous.png){: style="width:320px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/material-symmetries/icesym-ortho-viscous.png){: style="width:340px"} |

## Grain homogenization 

### Taylor&mdash;Sachs

Calculating $E_{ij}$ using (1) for a given CPO requires an *effective* rheology that takes the microstructure into account.

In the simplest case, polycrystals may be regarded as an ensemble of interactionless grains (monocrystals), subject to either a homogeneous stress field over the polycrystal scale:

$$
{\boldsymbol\tau}' = {\boldsymbol\tau}
,
\qquad\qquad \text{(Sachs's hypothesis)}
$$

or a homogeneous stain-rate field:

$$
{\dot{\boldsymbol\epsilon}}' = {\dot{\boldsymbol\epsilon}} 
,
\qquad\qquad \text{(Taylor's hypothesis)}
$$

where ${\boldsymbol\tau}'$ and ${\dot{\boldsymbol\epsilon}}'$ are the *microscopic* (grain-scale) stress and strain-rate tensors, respectively.

The effective rheology can then be approximated as the ensemble-averaged monocrystal rheology for either case:

$$
{\dot{\boldsymbol\epsilon}}^{\mathrm{Sachs}} = \langle {\dot{\boldsymbol\epsilon}}'({\boldsymbol\tau}') \rangle = \langle {\dot{\boldsymbol\epsilon}}'({\boldsymbol\tau}) \rangle
,
\qquad\qquad \text{(Sachs homogenization)}
$$

$$
\qquad
{\dot{\boldsymbol\epsilon}}^{\mathrm{Taylor}} = \langle {\boldsymbol\tau}'({\dot{\boldsymbol\epsilon}}') \rangle^{-1} = \langle {\boldsymbol\tau}'({\dot{\boldsymbol\epsilon}}) \rangle^{-1}
,
\qquad \text{(Taylor homogenization)}
$$

where $\langle \cdot \rangle^{-1}$ inverts the tensorial relationship.

If a linear combination of the two homogenizations is considered, equation (1) can be approximated as 

$$
E_{ij} = (1-\alpha) \frac{{\bf e}_i \cdot {\dot{\boldsymbol\epsilon}}^{\mathrm{Sachs}}({\boldsymbol\tau}) \cdot {\bf e}_j}
{{\bf e}_i \cdot {\dot{\boldsymbol\epsilon}}^{\mathrm{Sachs}}_{\mathrm{iso}}({\boldsymbol\tau}) \cdot {\bf e}_j }
+ {\alpha} \frac{ {\bf e}_i \cdot {\dot{\boldsymbol\epsilon}}^{\mathrm{Taylor}}({\boldsymbol\tau}) \cdot {\bf e}_j }
{ {\bf e}_i \cdot {\dot{\boldsymbol\epsilon}}^{\mathrm{Taylor}}_{\mathrm{iso}}({\boldsymbol\tau}) \cdot {\bf e}_j }
,
$$

or simply

$$
E_{ij} = (1-\alpha)E_{ij}^{\mathrm{Sachs}} + {\alpha}E_{ij}^{\mathrm{Taylor}} ,
$$

where $\alpha\in[0;1]$ is a free parameter.

!!! warning "Grain parameters"
    The grain viscous parameters used for homogenization should be understood as the *effective* values needed to reproduce deformation experiments on polycrystals; they are not the values derived from experiments on single crystals.

### Azuma&mdash;Placidi

🚧 *Not yet documented.*
