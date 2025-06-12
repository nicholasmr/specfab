# Representation

The Crystallographic Preferred Orientation (CPO) of polycrystalline materials is represented by the distribution of crystal axes in orientation space, $S^2$, weighted by grain size.

Supported grain symmetries, are:

| Grain symmetry | CPO components | Definition |
| --- | --- | --- | 
| Transversely isotropic | $n(\theta,\phi)$                  | Distribution of slip-plane normals |
| Orthotropic            | $n(\theta,\phi),\,b(\theta,\phi)$ | Distribution of slip-plane normals and slip directions |

Depending on which crystallographic slip system is preferentially activated, $n(\theta,\phi)$ and $b(\theta,\phi)$ refer to the distributions of different crystallographic axes.

!!! note "Glacier ice"

    | <center><div style="width:180px">Polycrystalline ice</div></center> | <center><div style="width:180px">Ensemble of slip elements</div></center> | |
    | :- | :- | :- |
    | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polyice-iso.png){: style="width:200px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/slip-plane/polycrystal-disk.png){: style="width:200px"} | <div style="width:200px;">Since ice grains are approximately<br> transversely isotropic, tracking<br> $n(\theta,\phi)$ (the $c$-axis distribution)<br> is sufficient for representing<br> the CPO of glacier ice. </div> | 

!!! tip "Olivine"

    | <center><div style="width:180px">Polycrystalline olivine</div></center> | <center><div style="width:180px">Ensemble of slip elements</div></center> | | 
    | :- | :- | :- |
    | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/polyoli-iso.png){: style="width:200px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/slip-plane/polycrystal-plane.png){: style="width:200px"} | <div style="width:200px;">For orthotropic grains such as<br> olivine, both $n(\theta,\phi)$ and $b(\theta,\phi)$<br> distributions must be tracked<br> to represent the CPO.<br><br> <u>*Notice*</u>: $n$ and $b$ represent the<br> distributions of a particular<br> crystallographic axes, depending<br> on fabric type (A&mdash;E type).</div> | 
   
## Background


In *specfab*, CPOs are  defined as the orientation distribution of grains without any reference to grain topology, grain shape, or the spatial arrangement of grains. 
More precisely, this means statistically describing the directions in which the slip-plane normal and slip direction of grains ($n$- and $b$-axes) are pointing, while also taking into account how massive each grain is (weighted by grain size). 


The *mass-density orientation distribution function* $\rho^{*}({\bf n},{\bf b})$ describes how grain mass is distributed in the space of possible grain orientations ([Faria, 2006](https://royalsocietypublishing.org/doi/abs/10.1098/rspa.2005.1610)), where ${\bf n}(\theta,\phi)$ and ${\bf b}(\theta,\phi)$ are arbitrary $n$- and $b$-axes (radial unit vectors). 
Since $\rho^{*}$ is by definition normalized by the polycrystal volume, integrating $\rho^{*}$ over all possible grain orientations recovers the mass density of the polycrystal (e.g., $\rho=917$ kg/m$^3$ for glacier ice):
$$ 
\rho = \int_{S^2} \int_{S^2} \rho^{*}({\bf n},{\bf b}) \,\mathrm{d}^2{\bf b}\, \mathrm{d}^2{\bf n} ,
$$ 
where integration is restricted to the surface of the unit sphere $S^2$ and $\mathrm{d}^2{\bf n}=\sin(\theta)\,\mathrm{d}{\theta}\,\mathrm{d}{\phi}$ is the infinitesimal solid angle in spherical coordinates (similarly for ${\bf b}$).

!!! warning "Ambiguity"

    Notice that defining the CPO in this way introduces some ambiguity. The contribution to $\rho^{*}({\bf n},{\bf b})$ from two grains with identical mass $m$ and orientation is indistinguishable from a single grain with the same orientation but twice the mass, $2m$. Nonetheless, how well $\rho^{*}({\bf n},{\bf b})$ represents a CPO should ultimately be judged by whether it contains the information necessary to calculate CPO-induced properties to some desired accuracy (say, bulk mechanical anisotropy); not by how simplified it is to disregard the spatial (topological) information of grains.

**Glacier ice**<br>
In the case of ice, *specfab* treats for simplicity all monocrystal properties as isotropic in the basal plane (transverse isotropy). 
This popular approach simplifies the problem significantly: since it does not matter in which direction $b$-axes (crystal $a$-axes) point, there is no need to track how they are distributed. 
Let therefore $n(\theta,\phi)$ be the corresponding normalized, *marginal* distribution function of the grain mass-density in the space of possible $n$-axis orientations (crystal $c$-axis orientations):
$$ 
n(\theta,\phi) = \frac{\int_{S^2} \rho^{*}({\bf n},{\bf b})\,\mathrm{d}^2{\bf b}}{\rho} .
$$

**Olivine**<br>
More complicated minerals like olivine are represented by also tracking the marginal distribution of slip directions:
$$ 
b(\theta,\phi) = \frac{\int_{S^2} \rho^{*}({\bf n},{\bf b})\,\mathrm{d}^2{\bf n}}{\rho} .
$$
If needed, the joint distribution function $\rho^{*}({\bf n},{\bf b})$ can be estimated from its marginal distributions following the identity for conditional probability density functions [and some assumptions](https://doi.org/10.1029/2024GC011831). 

## Mass or number density distributions?
The normalized, marginal distributions $n(\theta,\phi)$ and $b(\theta,\phi)$ are typically referred to as Mass Orientation Distribution Functions (MODFs) or Orientation Mass Distributions (OMDs).

In the literature, *number-density* distributions, rather than mass-density, frequently appear. 
In this case, $n(\theta,\phi)/N$ and $b(\theta,\phi)/N$ are referred to as the normalized Orientation Distribution Functions (ODFs) of slip-plane normals and slip directions, respectively, where $N = \int_{S^2} n(\theta,\phi) \,\mathrm{d}^2{\bf n}$ is the total number of grains. 

In the [models of crystal processes](fabdyn-matrix-model.md) available in *specfab*, the normalization is conserved and the two views give, in effect, the same result. 
The mass-density interpretation representation rests, however, on stronger physical grounds, since mass is conserved whereas grain numbers are not.

## Harmonic expansion

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/harmonic-expansion/harmonic-expansion.png#center){: style="width:630px"}

The distributions $n(\theta,\phi)$ and $b(\theta,\phi)$ are represented as an expansion series in spherical harmonics

$$ 
n(\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}n_{l}^{m}Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip-plane normals)},
$$

$$
b(\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}b_{l}^{m}Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip directions)}.
$$


The CPO state is therefore described by the *state vectors* of complex-valued expansion coefficients

$$
{\bf s}_n = [n_0^0,n_2^{-2},n_2^{-1},n_2^{0},n_2^{1},n_2^{2},n_4^{-4},\cdots,n_4^{4},\cdots,n_L^{-L},\cdots,n_L^{L}] \quad\text{($n$ state vector)}, 
$$

$$
{\bf s}_b = [b_0^0,b_2^{-2},b_2^{-1},b_2^{0},b_2^{1},b_2^{2},b_4^{-4},\cdots,b_4^{4},\cdots,b_L^{-L},\cdots,b_L^{L}] \quad\text{($b$ state vector)}, 
$$

where the magnitude and complex phase of each coefficient determines the size and rotation that a given harmonic contributes with.

### Reduced form

Not all expansion coefficients are independent for real-valued expansion series, but must fulfill

$$ 
n_l^{-m}=(-1)^m(n_l^m)^* ,
$$

$$ 
b_l^{-m}=(-1)^m(b_l^m)^* .
$$

This can be taken advantage of for large problems where many (e.g. gridded) CPOs must be stored in memory, thereby reducing the dimensionality of the problem. 
The vectors of reduced expansion coefficients are defined as

$$
\tilde{{\bf s}}_n = [n_0^0,n_2^{0},n_2^{1},n_2^{2},n_4^{0},\cdots,n_4^{4},\cdots,n_L^{0},\cdots,n_L^{L}] \quad\text{(reduced $n$ state vector)}.
$$

$$
\tilde{{\bf s}}_b = [b_0^0,b_2^{0},b_2^{1},b_2^{2},b_4^{0},\cdots,b_4^{4},\cdots,b_L^{0},\cdots,b_L^{L}] \quad\text{(reduced $b$ state vector)}.
$$

Converting between full and reduced forms is done as follows (similarly for $\tilde{{\bf s}}_b$):

```python
--8<-- "docs/snippets/reduced-form.py"
```

