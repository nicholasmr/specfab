# CPO representation

CPOs are represented by their distributions of crystallographic axes in orientation space ($S^2$). 
Supported grain symmetry groups are:

| Grain symmetry | CPO components | Definition |
| --- | --- | --- | 
| Transversely isotropic | $n(\theta,\phi)$                  | Distribution of slip-plane normals |
| Orthotropic            | $n(\theta,\phi),\,b(\theta,\phi)$ | Distribution of slip-plane normals and slip directions |

Thus, depending on which crystallographic slip system is preferentially activated, $n(\theta,\phi)$ and $b(\theta,\phi)$ may refer to the distributions of different crystallographic axes.

!!! note "Glacier ice"

    Since ice grains are approximately transversely isotropic, tracking $n(\theta,\phi)$ (the $c$-axis distribution) is sufficient for representing the CPO.

    | <center>Polycrystalline ice</center> | <center>Ensemble of slip elements</center> |
    | :- | :- |
    | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polycrystal-ice.png){: style="width:250px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/slip-plane/polycrystal-disk.png){: style="width:250px"} |

!!! note "Olivine"

    For orthotropic grains such as olivine, both $n(\theta,\phi)$ and $b(\theta,\phi)$ distributions must be tracked to represent the CPO.

    | <center>Polycrystalline olivine</center> | <center>Ensemble of slip elements</center> |
    | :- | :- |
    | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/polycrystal.png){: style="width:250px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/slip-plane/polycrystal-plane.png){: style="width:250px"} |

    Note that $n(\theta,\phi)$ and $b(\theta,\phi)$ represent the distributions of particular crystallographic axes (${\bf m}'_i$) depending on fabric type (A&mdash;E type).
    
   
## Definitions


CPOs of ice and olivine are here defined as the size and orientation distribution of grains, without any reference to grain topology, grain shape, or the spatial arrangement of grains. 
In the general-most case, this means statistically describing the directions in which grain slip-plane normal and slip direction ($n$- and $b$-axes) are pointing while also taking into account how massive each grain is. 


The *mass-density orientation distribution function* $\rho^{*}({\bf n},{\bf b})$ describes how grain mass is distributed in the space of possible grain orientations, where ${\bf n}$ and ${\bf b}$ are arbitrary $n$- and $b$-axes (unit vectors). Since $\rho^{*}$ is taken to be normalized by the polycrystal volume, integrating $\rho^{*}$ over all possible grain orientations recovers the mass density of the polycrystal (e.g. $\rho=917$ kg/m$^3$ for ice):
$$ 
\rho = \int_{S^2} \int_{S^2} \rho^{*}({\bf n},{\bf b}) \,\mathrm{d}^2{\bf b}\, \mathrm{d}^2{\bf n} ,
$$ 
where integration is restricted to the surface of the unit sphere $S^2$, and $\mathrm{d}^2{\bf n}=\sin(\theta)\,\mathrm{d}{\theta}\,\mathrm{d}{\phi}$ is the infinitesimal solid angle in spherical coordinates (defined similarly for ${\bf a}$).

!!! warning "Ambiguity"

    Notice that defining CPOs in this way introduces some ambiguity. The contribution to $\rho^{*}({\bf n},{\bf b})$ from two grains with identical mass $m$ and orientation is indistinguishable from a single grain with the same orientation but twice the mass, $2m$. Nonetheless, how well $\rho^{*}({\bf n},{\bf b})$ can represent CPOs should ultimately be judged by whether it contains the information necessary to calculate CPO-induced properties to some desired accuracy (say, bulk mechanical anisotropy); not by how simplified it is to disregard spatial (topological) grain information etc. 

**Glacier ice**<br>
In the case of ice, specfab treats for simplicity approximate all monocrystal properties as isotropic in the basal plane (transverse isotropy). 
This popular approach simplifies the problem significantly: since it does not matter in which directions $b$-axes (crystal $a$-axes) point, it is not needed to keep track of how they are distributed.
Let therefore $n(\theta,\phi)$ be the corresponding distribution function of grain mass density in the space of possible $n$-axis (crystal $c$-axis) orientations:
$$ 
n(\theta,\phi) = \frac{\int_{S^2} \rho^{*}({\bf n},{\bf b})\,\mathrm{d}^2{\bf b}}{\rho} ,
$$
which is nothing but the normalized marginal distribution of ${\bf n}$.

**Olivine**<br>
Specfab represents more complicated minerals like olivine by, in addition, also tracking the marginal distribution of slip directions:
$$ 
b(\theta,\phi) = \frac{\int_{S^2} \rho^{*}({\bf n},{\bf b})\,\mathrm{d}^2{\bf n}}{\rho} .
$$
When needed, the joint distribution function $\rho^{*}({\bf n},{\bf b})$ is the back-constructed from its marginal distributions following the identity for conditional probability density functions [and some assumptions](https://doi.org/10.1029/2024GC011831). 

<!--

The orientation distribution function (ODF) of a given slip-system axis (crystallographic axis) $f\in \lbrace n,b\rbrace$ is defined as the normalized distribution

$$ 
\mathrm{ODF} = \frac{f(\theta,\phi)}{N} \quad\text{where}\quad N=\int_{S^2} f(\theta,\phi) \,\mathrm{d}\Omega .
$$

## Normalization

$n(\theta,\phi)$ may be understood either as the number density of grains with a given slip-plane normal orientation, or as the mass density fraction ([Faria, 2006](https://royalsocietypublishing.org/doi/abs/10.1098/rspa.2005.1610); [Richards et al., 2021](https://www.sciencedirect.com/science/article/abs/pii/S0012821X20306622)) of grains with a given slip-plane normal orientation; $\varrho^*(\theta,\phi)$ in literature.
The same goes for $b(\theta,\phi)$.

From specfab's point-of-view, the difference is a matter of normalization: since the models of [CPO evolution](cpo-dynamics-tranisotropic.md) (lattice rotation, DDRX, CDRX) conserve the normalization, the two views are effectively the same, not least because CPO-derived quantities depend on the normalized distributions (which are identical).
The mass-density-fraction interpretation rests, however, on stronger physical grounds as mass is conserved but grain numbers are not.

-->


## Harmonic expansion

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/harmonic-expansion/harmonic-expansion.png#center){: style="width:750px"}

The distributions $n(\theta,\phi)$ and $b(\theta,\phi)$ are represented as spherical harmonic expansion series:

$$ 
n(\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}n_{l}^{m}Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip-plane normals)},
$$

$$
b(\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}b_{l}^{m}Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip directions)}.
$$


The CPO state is thus described by the *state vectors* of complex-valued expansion coefficients

$$
{\bf s}_n = [n_0^0,n_2^{-2},n_2^{-1},n_2^{0},n_2^{1},n_2^{2},n_4^{-4},\cdots,n_4^{4},\cdots,n_L^{-L},\cdots,n_L^{L}] \quad\text{($n$ state vector)}, 
$$

$$
{\bf s}_b = [b_0^0,b_2^{-2},b_2^{-1},b_2^{0},b_2^{1},b_2^{2},b_4^{-4},\cdots,b_4^{4},\cdots,b_L^{-L},\cdots,b_L^{L}] \quad\text{($b$ state vector)}, 
$$

where the magnitude and complex phase of the coefficients determine the size and rotation of the contribution from the associated harmonic mode.

### Reduced form

Not all expansion coefficients are independent for real-valued expansion series, but must fulfill (likewise for $b$)

$$ 
n_l^{-m}=(-1)^m(n_l^m)^* .
$$

This can be taken advantage of for large problems where many (e.g. gridded) CPOs must be stored in memory, thereby effectively reducing the size of the problem. 
The vector of reduced expansion coefficients is defined as

$$
\tilde{{\bf s}}= [n_0^0,n_2^{0},n_2^{1},n_2^{2},n_4^{0},\cdots,n_4^{4},\cdots,n_L^{0},\cdots,n_L^{L}] \quad\text{(reduced state vector)}.
$$

Converting between full and reduced forms is done as follows:

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(2) # L=2 truncation is sufficient in this case

### Construct an arbitrary fabric
a2 = np.diag([0.1,0.2,0.7]) # arbitrary second-order structure tensor
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # determine l<=2 expansion coefficients of ODF
print('original:', nlm)

### Get reduced form of coefficient array, rnlm
rnlm_len = sf.get_rnlm_len() 
rnlm = np.zeros((rnlm_len), dtype=np.complex64) # array of reduced expansion coefficients
rnlm[:] = sf.nlm_to_rnlm(nlm, rnlm_len) # reduced form
print('reduced:', rnlm)

### Recover full form (nlm) from reduced form (rnlm)
nlm[:] = sf.rnlm_to_nlm(rnlm, nlm_len)
print('recovered:', nlm)
```

