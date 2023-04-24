# CPO representation

## Definition

CPOs are represented by the distribution(s) of crystallographic axes in orientation space, $S^2$.
<br>
Supported grain symmetry groups for modelling [CPO evolution](cpo-dynamics-tranisotropic.md) are

| Grain symmetry | CPO components | Interpretation |
| --- | --- | --- | 
| Transversely isotropic | $n(\theta,\phi)$                | Distribution of slip-plane normals |
| Orthotropic            | $n(\theta,\phi),b(\theta,\phi)$ | Distribution of slip-plane normals and slip directions |

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/slipplane.png){: style="width:180px"} 

!!! note
    Grain sizes or shapes are not modelled by specfab.

### Example

| <center>Polycrystalline ice</center> | <center>Polycrystalline olivine</center> |
| :- | :- |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/polycrystal.png){: style="width:220px"} | ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/polycrystal.png){: style="width:220px"} |
| $n(\theta,\phi)$ is the ${\bf c}$-axis distribution for <br>polycrystalline ice. | $n(\theta,\phi)$ and $b(\theta,\phi)$ are the distributions of particular <br>crystallographic axes (${\bf r}_i$) for polycrystalline olivine<br> depending on fabric type (A&mdash;E type). |

## Series expansion

CPOs are represented by expanding their distributions of crystallographic axes in terms of spherical harmonic expansion series.

E.g. for transversely isotropic grains where only $n(\theta,\phi)$ is relevant:
$$ 
n(\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}n_{l}^{m}Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip-plane normals)}.
$$

The orientation distribution function (ODF) is defined as the normalized distribution 

$$ 
\mathrm{ODF} = \frac{n(\theta,\phi)}{N} \quad\text{where}\quad N=\int_{S^2}{n} d\Omega=\sqrt{4\pi}n_0^0 .
$$

The array of complex-valued expansion coefficients, defining the CPO state, is 

$$
{\bf s} = [n_0^0,n_2^{-2},n_2^{-1},n_2^{0},n_2^{1},n_2^{2},n_4^{-4},\cdots,n_4^{4},\cdots,n_L^{-L},\cdots,n_L^{L}] \quad\text{(state vector)}.
$$

$$ $$ <!-- half space -->

!!! warning "Normalization"

    $n(\theta,\phi)$ may be understood either as the number density of grains with a given slip-plane normal orientation, or as the mass density fraction ([Faria, 2006](https://royalsocietypublishing.org/doi/abs/10.1098/rspa.2005.1610); [Richards et al., 2021](https://www.sciencedirect.com/science/article/abs/pii/S0012821X20306622)) of grains with a given slip-plane normal orientation.

    From specfab's point-of-view, the difference is a matter of normalization: since the models of [CPO evolution](cpo-dynamics-tranisotropic.md) (lattice rotation, DDRX, CDRX) conserve the normalization, the two views are effectively the same, not least because CPO-derived quantities depend on the normalized distributions (which are identical).
    The mass-density-fraction interpretation rests, however, on stronger physical grounds as mass is conserved but grain numbers are not.

## Structure tensors

The 2nd-, 4th- and 6th-order structure tensors (`a2`, `a4`, `a6`) are the average outer products of a crystallographic axis ${\bf r}$ (here assuming equal grain weight/size/mass for simplicity)

$$ 
{\bf a}^{(2)}=\frac{1}{N}\sum_i^N {\bf r}_i\otimes{\bf r}_i 
,\quad 
{\bf a}^{(4)}=\frac{1}{N}\sum_i^N ({\bf r}_i\otimes)^4
,\quad 
{\bf a}^{(6)}=\frac{1}{N}\sum_i^N ({\bf r}_i\otimes)^6
.
$$


Converting between spectral and tensorial representations is a linear problem in the sense that 

$$ {\bf a}^{(2)} = {\bf f}(\hat{n}_2^{m})
,\quad
{\bf a}^{(4)} = {\bf g}(\hat{n}_0^0, \hat{n}_2^{m}, \hat{n}_4^{m}) 
,\quad
{\bf a}^{(6)} = {\bf h}(\hat{n}_2^{m}, \hat{n}_4^{m}, \hat{n}_6^{m}) 
,
\qquad\text{(for all $m$)}
$$

where $\hat{n}_l^m = n_l^m/n_0^0$, and ${\bf f}$, ${\bf g}$, ${\bf h}$ are linear in their arguments.

### Example

The following example shows how to convert between the representations:

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(6) # L=6 truncation is sufficient in this case
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients

### a2 to nlm
a2 = np.diag([0.0,0.25,0.75]) # any second-order structure tensor (not necessarily diagonal)
nlm[:6] = sf.a2_to_nlm(a2) # Determine l=2 expansion coefficients of ODF (a2 is normalized)
a2 = sf.a2(nlm) # nlm back to a2
print('a2 is: ', a2)

### a4 to nlm
p = np.array([0,0,1]) # unidirectional fabric where all c-axes are aligned with the z-direction
a4 = np.einsum('i,j,k,l', p,p,p,p) # a4 for ODF = deltafunc(r-p) 
nlm[:15] = sf.a4_to_nlm(a4) # Determine l=2,4 expansion coefficients of ODF (a4 is normalized)
a4 = sf.a4(nlm) # nlm back to a4 
print('a4 is: ', a4)

### a6 to nlm
p = np.array([0,0,1]) # unidirectional fabric where all c-axes are aligned with the z-direction
a6 = np.einsum('i,j,k,l,m,n', p,p,p,p,p,p) # a6 for ODF = deltafunc(r-p) 
nlm[:] = sf.a6_to_nlm(a6) # Determine l=2,4,6 expansion coefficients of ODF (a6 is normalized)
a6 = sf.a6(nlm) # nlm back to a6
print('a6 is: ', a6)
```

## Constructing CPOs from measurements

The spectral expansion coefficients of any CPO may be determined from discrete measurements of crystallographic axes.
This requires constructing the corresponding structure tensors (for each crystallographic axis), from which the expansion coefficients may be derived.

###Example
```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(6) 

# Replace with your own array/list of measured c-axes
# caxes = [[c1x,c1y,c1z], [c2x,c2y,c2z], ...] 

# Determine sixth-order structure tensor, a6
a6 = np.zeros((3,3,3,3,3,3))
for c in caxes
    a6 += np.einsum('i,j,k,l,m,n', c,c,c,c,c,c) # sixth outer product of c with itself
a6 /= len(caxes) # normalize by number of c-axes (grains)

# Determine spectral expansion coefficients
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[:] = sf.a6_to_nlm(a6) # Determine l=2,4,6 expansion coefficients of ODF (a6 is normalized)
```

!!! note 
    Constructing `a6` is to be preferred over `a4` and `a2` since it contains more information on the fine-scale structure of the distribution; 
    that is, $l\leq 6$ expansion coefficients as opposed to $l\leq 4$ and $l\leq 2$ coefficients, respectively.


## Reduced form

Not all expansion coefficients are independent for real-valued expansion series, but must fulfill

$$ 
n_l^{-m}=(-1)^m(n_l^m)^* .
$$

This can be taken advantage of for large problems where many (e.g. gridded) CPOs need to be stored in memory, thereby effectively reducing the size of the problem. 
The array of reduced expansion coefficients is defined as

$\qquad$ `rnlm` $= [n_0^0,n_2^{0},n_2^{1},n_2^{2},n_4^{0},\cdots,n_4^{4},\cdots,n_L^{0},\cdots,n_L^{L}] \quad\text{(reduced state vector)}$

### Example 

Converting between full and reduced forms:

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(2) # L=2 truncation is sufficient in this case

# Construct an arbitrary fabric
a2 = np.diag([0.1,0.2,0.7]) # any second-order structure tensor (not necessarily diagonal)
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[0:6] = sf.a2_to_nlm(a2) # l=2 expansion coefficients for corresponding ODF (a2 is normalized)
print('original:', nlm)

# Get reduced form of coefficient array, rnlm
rnlm_len = sf.get_rnlm_len() 
rnlm = np.zeros((rnlm_len), dtype=np.complex64) # array of reduced expansion coefficients
rnlm[:] = sf.nlm_to_rnlm(nlm, rnlm_len) # reduced form
print('reduced:', rnlm)

# Recover full form (nlm) from reduced form (rnlm)
nlm[:] = sf.rnlm_to_nlm(rnlm, nlm_len)
print('recovered:', nlm)
```

