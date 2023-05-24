# Structure tensors

The 2nd-, 4th- and 6th-order structure tensors (`a2`, `a4`, `a6`) are the average outer products of a crystallographic axis ${\bf r}_i$ (here assuming equal grain weight/size/mass for simplicity)

$$ 
{\bf a}^{(2)}=\frac{1}{N}\sum_i^N {\bf r}_i\otimes{\bf r}_i 
,\quad 
{\bf a}^{(4)}=\frac{1}{N}\sum_i^N ({\bf r}_i\otimes)^4
,\quad 
{\bf a}^{(6)}=\frac{1}{N}\sum_i^N ({\bf r}_i\otimes)^6
,
$$

where $N$ is the total number of grains. 
In terms of the spectral expansion $n(\theta,\phi)$, the structure tensors are

$$ 
{\bf a}^{(k)}=\frac{1}{N} \int_{S^2} ({\bf r}_i\otimes)^k n(\theta,\phi) \, \mathrm{d}\Omega
,
$$

where $N=\int_{S^2} n(\theta,\phi) \, \mathrm{d}\Omega$ and $\mathrm{d}\Omega = \mathrm{d}\theta \mathrm{d}\phi$ is the infinitesimal solid angle.


Converting between spectral and tensorial representations is a linear problem in the sense that 

$$ {\bf a}^{(2)} = {\bf f}(\hat{n}_2^{m})
,\quad
{\bf a}^{(4)} = {\bf g}(\hat{n}_0^0, \hat{n}_2^{m}, \hat{n}_4^{m}) 
,\quad
{\bf a}^{(6)} = {\bf h}(\hat{n}_2^{m}, \hat{n}_4^{m}, \hat{n}_6^{m}) 
,
\qquad\text{(for all $m$)}
$$

where ${\bf f}$, ${\bf g}$, ${\bf h}$ are linear in their arguments, and 

$$
\hat{n}_l^m \equiv n_l^m/n_0^0
.
$$

## Example

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

# Constructing CPOs from measurements

The spectral expansion coefficients of any CPO may be determined from discrete measurements of crystallographic axes.
This requires constructing the corresponding structure tensors (for each crystallographic axis), from which the expansion coefficients may be derived.

## Example
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

