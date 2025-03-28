# Radar-derived physical properties of glacier ice  

## Introduction

The dielectric permittivity tensor of a single ice crystal is approximately transversely isotropic w.r.t. the crystal $c$-axis:
$$
\epsilon_{ij}' = (2\epsilon_{\perp}' + \epsilon_{\parallel}') \frac{\delta_{ij}}{3}
+ (\epsilon_{\parallel}'-\epsilon_{\perp}') \left(c_i c_j - \frac{\delta_{ij}}{3} \right),
$$
where $\epsilon_{\parallel}'$ and $\epsilon_{\perp}'$ are the components parallel and perpendicular to the $c$-axis, respectively, which depend on ice temperature and EM-wave frequency ([Fujita et al., 2000](https://eprints.lib.hokudai.ac.jp/dspace/bitstream/2115/32469/1/P185-212.pdf)).

For wave lengths much longer than the average grain size, the bulk permittivity tensor of polycrystalline ice may be approximated as the grain-average permittivity tensor, $\epsilon_{ij} \simeq \langle \epsilon_{ij}' \rangle$, constructed by averaging over all grain orientations (over the CPO) assuming grain sizes are uncorrelated with orientation:

$$
\langle \epsilon_{ij}' \rangle = 
(2\epsilon_{\perp}' + \epsilon_{\parallel}') \frac{\delta_{ij}}{3}
+ (\epsilon_{\parallel}'-\epsilon_{\perp}') \left(\langle c_i c_j \rangle - \frac{\delta_{ij}}{3} \right)
,
$$

where $\langle c_i c_j \rangle$ is the second-order [structure tensor](cpo-structuretensors.md), defined as

$$ 
\langle c_i c_j \rangle = 
\frac{1}{N}\sum_{k=1}^N { c_i^{(k)} c_j^{(k)} }.
$$

Thus, because the bulk permittivity tensor $\epsilon_{ij}$ can be inferred from [EM-wave speeds](wavepropagation-electromagnetic.md) and radar return-power anomalies, so can $\langle c_i c_j \rangle$.

## Radar measurements $\rightarrow$ CPO

A useful approximation over large parts of ice sheets is that $\langle c_i c_j \rangle$ has a vertical eigenvector, in which case the Cartesian components are

$$
\langle c_i c_j \rangle 
=
\left[\begin{matrix}
a_{xx} & a_{xy} & 0\\ 
a_{xy} & a_{yy} & 0\\ 
0 & 0 & a_{zz}
\end{matrix}\right]
.
$$

Let us consider the usual case where the difference in horizontal eigenvalues of $\langle c_i c_j \rangle$,

$$
\Delta \lambda = \lambda_2 - \lambda_1,
$$

can be inferred from ice-penetrating radar, where ${\bf m}_1$ and ${\bf m}_2$ are the corresponding horizontal eigenvectors and eigenvalues are sorted such that $\lambda_1 \leq \lambda_2$.
It follows that the structure tensor, posed in its eigenframe (${\bf m}_1, {\bf m}_2, {\bf z}$), is

$$
\langle c_i c_j \rangle 
=
\left[\begin{matrix}
\lambda_1  & 0 & 0 \\ 
0 & \lambda_1 + \Delta\lambda  & 0 \\ 
0 & 0 & 1 - \Delta \lambda - 2\lambda_1
\end{matrix}\right]
,
$$

where the identity $\operatorname{tr}(\langle c_i c_j \rangle) = 1$ was used.

### Gerber's approximation 

Since $\lambda_1$ is unknown, the problem can be closed by making different assumptions about $\lambda_1$ given the local/upstream flow regime, such as proposed by [Gerber et al. (2023)](https://www.nature.com/articles/s41467-023-38139-8).

Suppose $\Delta\lambda$ is measured in region where $c$-axes are, to a good approximation, suspected to be distributed on the ${\bf m}_2$&mdash;${\bf z}$ plane because the smallest eigenvalue is vanishing, $\lambda_1 \rightarrow 0$.
In this case, $\Delta \lambda = 0$ represents a perfect single-maximum along ${\bf z}$, $\Delta \lambda = 0.5$ a perfect girdle in the ${\bf m}_2$&mdash;${\bf z}$ plane, and $\Delta \lambda = 1$ a perfect single-maximum along ${\bf m}_2$, respectively:

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/radar-PP-figs/plane-CPOs.png){: style="width:650px"}

## CPO $\rightarrow$ Enhancement factors

If $\langle c_i c_j \rangle$ can be inferred from radar sounding following the above method, so can the bulk strain-rate enhancement factors, $E_{ij}$, in the same eigenframe (i.e. [eigenenhancements](enhancements-strainrate.md)).

The eigenenhancements depend, however, also on the fourth-order structure tensor, $\langle c_i c_j c_k c_l \rangle$, but the bulk permittivity $\epsilon_{ij}$ is insensitive to $\langle c_i c_j c_k c_l \rangle$.
To overcome this, a simple empirical correlation is adopted that allows determining $\langle c_i c_j c_k c_l \rangle$ given $\langle c_i c_j\rangle$ if the CPO is approximately rotationally symmetric.

### Correlation between $\langle c_i c_j c_k c_l \rangle$ and $\langle c_i c_j\rangle$

If the CPO symmetry axis is rotated into the vertical direction, $\langle c_i c_j\rangle$ depends only on the normalized [spectral component](cpo-representation.md) $\hat{n}_2^0 = n_2^0/n_0^0:$

$$
\langle c_i c_j\rangle = \frac{\delta_{ij}}{3} +  \frac{2\sqrt{5}}{15} \hat{n}_2^0
\left[\begin{matrix}
-1/2 & 0 & 0 \\ 
0  & -1/2  & 0 \\ 
0 & 0 & 1
\end{matrix}\right]
,
$$

and $\langle c_i c_j c_k c_l \rangle$ only on $\hat{n}_2^0$ and $\hat{n}_4^0 = n_4^0/n_0^0$ (not shown).
The figure below shows the empirical correlation between these two components based on ice-core samples.

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/research/state-space/ice/state-space-empcorr.png){: style="width:570px"}

Thus, if $\hat{n}_2^0$ is extracted from $\langle c_i c_j\rangle$ in this frame, $\hat{n}_4^0$ can be derived and hence $\langle c_i c_j c_k c_l \rangle$ constructed.
To pose the CPO in the original, unrotated eigenframe (${\bf m}_1, {\bf m}_2, {\bf z}$), the resulting expansion series is finally rotated back, allowing eigenenhancements to easily be calculated using specfab.

### Code example

The following code demonstrates how to take each step with specfab:

```python
import numpy as np
from scipy.spatial.transform import Rotation
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) # L=4 is sufficient here

### Determine <c_i c_j> from radar-derived Delta lambda
l1 = 0   # lambda_1 = 0 (Gerber's approximation)
dl = 0.5 # Delta lambda = lambda_2 - lambda_1
a2 = np.diag([l1, l1+dl, 1-dl-2*l1]) # second-order structure tensor, <c_i c_j>, in eigenframe
m1, m2, z = np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1]) # eigenvectors

### Rotate <c_i c_j> into a rotationally-symmetric frame about z
Rm1 = Rotation.from_rotvec(np.pi/2 * m1).as_matrix() # Rotate 90 deg about m1 eigenvector
Rm2 = Rotation.from_rotvec(np.pi/2 * m2).as_matrix() # Rotate 90 deg about m2 eigenvector
if dl < 0.4:         a2_vs = a2 # Already in rotationally-symmetric frame about z
if 0.4 <= dl <= 0.6: a2_vs = np.matmul(Rm2,np.matmul(a2,Rm2.T)) # Rotate vertical (m2--z) girdle into horizontal (m1--m2) girdle
if dl > 0.6:         a2_vs = np.matmul(Rm1,np.matmul(a2,Rm1.T)) # Rotate horizontal (m2) single-maximum into vertical (z) single-maximum

### Determine \hat{n}_4^0 (= n_4^0/n_0^0) from \hat{n}_2^0 (= n_2^0/n_0^0) in rotationally-symmetric frame about z
nhat20 = (a2_vs[2,2]- 1/3)/(2/15*np.sqrt(5)) # azz -> nhat20
nhat40 = sf.nhat40_empcorr_ice(nhat20)[0]

### Construct nlm (spectral CPO state vector) in rotationally-symmetric frame about z
nlm_vs = np.zeros(nlm_len, dtype=np.complex128) 
n00 = 1/np.sqrt(4*np.pi) # only grain-number normalized distribution is known, so must integrate to 1 over S^2.
nlm_vs[0]  = n00
nlm_vs[3]  = nhat20*n00
nlm_vs[10] = nhat40*n00

### Rotate spectral CPO state back to origional (m1,m2,z) eigenframe 
if dl < 0.4:         nlm = nlm_vs # Already in vertical symmetric frame
if 0.4 <= dl <= 0.6: nlm = sf.rotate_nlm(nlm_vs, -np.pi/2, 0) # Rotate horizontal (m1--m2) girdle back into vertical (m2--z) girdle
if dl > 0.6:         nlm = sf.rotate_nlm(sf.rotate_nlm(nlm_vs, -np.pi/2, 0), 0 ,-np.pi/2) # Rotate vertical (z) single-maximum back into horizontal (m2) single-maximum

### Calculate eigenenhancements
# Transversely isotropic monocrystal parameters for ice (Rathmann & Lilien, 2021)
n_grain   = 1        # Power-law exponent: n=1 => linear grain rheology, nonlinear (n>1) is unsupported
Eij_grain = (1, 1e3) # Grain eigenenhancements (Ecc,Eca) for compression along c-axis (Ecc) and for shear parallel to basal plane (Eca)
alpha     = 0.0125   # Taylor--Sachs weight
# Tuple of eigenenhancements (bulk enhancement factors w.r.t. m1, m2, z)
e1, e2, e3 = m1, m2, z
Eij = sf.Eij_tranisotropic(nlm, e1,e2,e3, Eij_grain,alpha,n_grain) # (E_{m1,m1},E_{m2,m2},E_{zz},E_{m2,z),E_{m1,z},E_{m1,m2})
# To calculate bulk enhancement factors w.r.t. other axes of deformation/stress, change (e1,e2,e3) accordingly.
```

For reference, the below plots show the different CPOs at each step for $\Delta\lambda=0.5$ and $\Delta\lambda=1$.

!!! info "$\Delta\lambda = 0.5$"

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/radar-PP-figs/code-example-output-dl0.5.png){: style="width:650px"}

!!! info "$\Delta\lambda = 1.0$"

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/radar-PP-figs/code-example-output-dl1.0.png){: style="width:650px"}


