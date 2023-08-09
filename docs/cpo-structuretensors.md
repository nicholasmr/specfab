# Structure tensors

The $k$-th order structure tensor is the average $k$-th repeated outer product of a crystallographic axis with itself (vector moments).
For example, in the case of a discrete ensemble of ${\bf c}$ axes they are

$$ 
{\bf a}^{(k)}=\frac{1}{N}\sum_i^N ({\bf c}_i\otimes)^k,
$$

where $N$ is the total number of grains, assuming equal grain weight (i.e. mass) for simplicity.
Alternatively, if the distribution function of ${\bf c}$ axes is known, $n(\theta,\phi)$, the structure tensors are

$$ 
{\bf a}^{(k)}=\frac{1}{N} \int_{S^2} (\hat{{\bf r}}\otimes)^k n(\theta,\phi) \, \mathrm{d}\Omega
,
$$

where $\mathrm{d}\Omega = \sin(\theta) \mathrm{d}\theta \mathrm{d}\phi$ is the infinitesimal solid angle, $\hat{{\bf r}}(\theta,\phi)$ is the radial unit vector, and $N=\int_{S^2} n(\theta,\phi) \, \mathrm{d}\Omega$.

## Convert to spectral 

Converting between spectral and tensorial representations is a linear problem in the sense that 

$$
{\bf a}^{(k)} = {\bf f}(\hat{n}_2^{m}, \hat{n}_4^{m}, \cdots, \hat{n}_k^{m}) 
,
\qquad\text{(for all $m$)}
$$

where ${\bf f}$ is linear in its arguments, and 

$$
\hat{n}_l^m = n_l^m/n_0^0
.
$$

In the case of ${\bf a}^{(2)}$ the relation is simple:

$$
{\bf a}^{(2)} = \frac{{\bf I}}{3} + \sqrt{\frac{2}{15}}
\left[\begin{matrix}
\operatorname{Re}[\hat{n}_2^2] - \dfrac{1}{2}\sqrt{\dfrac{2}{3}} \hat{n}_2^0 & -\operatorname{Im}[\hat{n}_2^2] & -\operatorname{Re}[\hat{n}_2^1] \\ 
 & -\operatorname{Re}[\hat{n}_2^2] - \dfrac{1}{2}\sqrt{\dfrac{2}{3}} \hat{n}_2^0  & \operatorname{Im}[\hat{n}_2^1] \\ 
\mathrm{sym.} &  & \sqrt{\dfrac{2}{3}} \hat{n}_2^0
\end{matrix}\right]
,
$$

but for higher-order structure tensors the expressions are long (not shown).

The following code example shows how to convert between the representations:

```python
import numpy as np
from specfabpy import specfab as sf
L = 8
lm, nlm_len = sf.init(L)
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients

### a2 to nlm
a2 = np.diag([0.0,0.25,0.75]) # arbitrary second-order structure tensor
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # determine l<=2 expansion coefficients of ODF
a2 = sf.a2(nlm) # nlm back to a2
print('a2 is: ', a2)

### a4 to nlm
p = np.array([0,0,1]) # unidirectional CPO
a4 = np.einsum('i,j,k,l', p,p,p,p) # a4 for ODF = deltafunc(r-p) 
nlm[:sf.L4len] = sf.a4_to_nlm(a4) # determine l<=4 expansion coefficients of ODF
a4 = sf.a4(nlm) # nlm back to a4 
print('a4 is: ', a4)

### a6 to nlm
p = np.array([0,0,1]) # unidirectional CPO
a6 = np.einsum('i,j,k,l,m,n', p,p,p,p,p,p) # a6 for ODF = deltafunc(r-p) 
nlm[:sf.L6len] = sf.a6_to_nlm(a6) # determine l<=6 expansion coefficients of ODF
a6 = sf.a6(nlm) # nlm back to a6
print('a6 is: ', a6)
```

## Construct from measurements

The spectral expansion coefficients of any CPO may be determined from discrete measurements of crystallographic axes.
This requires constructing the corresponding structure tensors (for each crystallographic axis), from which the expansion coefficients may be derived.

```python
import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(8) 

### Replace with your own array/list of measured c-axes
# caxes = [[c1x,c1y,c1z], [c2x,c2y,c2z], ...] 

### Determine sixth-order structure tensor, a6
a6 = np.zeros((3,3,3,3,3,3))
for c in caxes
    a6 += np.einsum('i,j,k,l,m,n', c,c,c,c,c,c) # sixth outer product of c-axis with itself
a6 /= len(caxes) # normalize by number of c-axes (grains)

### Determine spectral expansion coefficients
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[:sf.L6len] = sf.a6_to_nlm(a6) # determine l<=6 expansion coefficients of ODF 
```

Note that constructing `a6` is to be preferred over `a4` and `a2` since it contains more information on the fine-scale structure of the distribution; 
that is, $l\leq 6$ expansion coefficients as opposed to $l\leq 4$ and $l\leq 2$ coefficients, respectively.

