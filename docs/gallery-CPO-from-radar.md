# Inferring CPO from radar

The dielectric permittivity tensor of a single ice crystal is often taken to be transversely isotropic w.r.t. the crystal $c$-axis:

$$
\boldsymbol{\epsilon}' = (2\epsilon_{\perp}' + \epsilon_{\parallel}') \frac{{\bf I}}{3}
+ (\epsilon_{\parallel}'-\epsilon_{\perp}') \left({\bf c}^2 - \frac{{\bf I}}{3} \right),
$$

where $\epsilon_{\parallel}'$ and $\epsilon_{\perp}'$ are the components parallel and perpendicular to the $c$-axis, respectively, which depend on ice temperature and EM wave frequency ([Fujita et al., 2000](https://eprints.lib.hokudai.ac.jp/dspace/bitstream/2115/32469/1/P185-212.pdf)).

For wave lengths much longer than the average grain size, the bulk permittivity tensor $\boldsymbol{\epsilon}$ of *polycrystalline* ice may be approximated as the grain-average permittivity tensor, constructed by averaging $\boldsymbol{\epsilon}'$ over all grain orientations (over the CPO) 

$$
\boldsymbol{\epsilon} \simeq \langle \boldsymbol{\epsilon}' \rangle = 
(2\epsilon_{\perp}' + \epsilon_{\parallel}') \frac{{\bf I}}{3}
+ (\epsilon_{\parallel}'-\epsilon_{\perp}') \left(\langle {\bf c}^2 \rangle - \frac{{\bf I}}{3} \right)
,
$$

where $\langle {\bf c}^2\rangle$ is the second-order [structure tensor](cpo-structuretensors.md), identical to ${\bf a}^{(2)}$.

Thus, because the bulk permittivity tensor $\boldsymbol{\epsilon}$ can be inferred from [EM wave speeds](waveprop-electromagnetic.md) or radar return-power anomalies, so can $\langle {\bf c}^2 \rangle$.

## Radar measurements $\rightarrow$ CPO

A useful approximation over large parts of ice sheets is that $\langle {\bf c}^2 \rangle$ has a vertical eigenvector, in which case the Cartesian components are

$$
\langle {\bf c}^2 \rangle =
\left[\begin{matrix}
a_{xx} & a_{xy} & 0\\ 
a_{xy} & a_{yy} & 0\\ 
0 & 0 & a_{zz}
\end{matrix}\right]
.
$$

Consider then the case where the difference in horizontal eigenvalues of $\langle {\bf c}^2 \rangle$,

$$
\Delta \lambda = \lambda_2 - \lambda_1,
$$

can be [inferred from ice-penetrating radar](waveprop-electromagnetic.md), where ${\bf m}_1$ and ${\bf m}_2$ are the corresponding horizontal eigenvectors and eigenvalues are sorted such that $\lambda_1 \leq \lambda_2$.
It follows that the structure tensor, posed in its eigenframe (${\bf m}_1, {\bf m}_2, {\bf z}$), is

$$
\langle {\bf c}^2 \rangle
=
\left[\begin{matrix}
\lambda_1  & 0 & 0 \\ 
0 & \lambda_1 + \Delta\lambda  & 0 \\ 
0 & 0 & 1 - \Delta \lambda - 2\lambda_1
\end{matrix}\right]
,
$$

where the identity $\operatorname{tr}(\langle {\bf c}^2 \rangle) = 1$ was used.

### 💡 Gerber's approximation 

Since $\lambda_1$ is unknown, the problem can be closed by making different assumptions about $\lambda_1$ given the local/upstream flow regime, such as proposed by [Gerber et al. (2023)](https://www.nature.com/articles/s41467-023-38139-8).

Suppose $\Delta\lambda$ is measured in region where $c$-axes are, to a good approximation, suspected to be distributed on the ${\bf m}_2$&mdash;${\bf z}$ plane because the smallest eigenvalue is vanishing, $\lambda_1 \rightarrow 0$.
In this case, $\Delta \lambda = 0$ represents a perfect single-maximum along ${\bf z}$, $\Delta \lambda = 0.5$ a perfect girdle in the ${\bf m}_2$&mdash;${\bf z}$ plane, and $\Delta \lambda = 1$ a perfect single-maximum along ${\bf m}_2$, respectively:

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/radar-PP-figs/plane-CPOs.png){: style="width:600px"}

## CPO $\rightarrow$ Enhancement factors

If $\langle {\bf c}^2 \rangle$ can be inferred from radar measurements following the above method, so can the bulk strain-rate [enhancement factors](enhancements-strainrate.md) $E_{ij}$ in the same frame.
The enhancements depend, however, also on the fourth-order structure tensor, $\langle {\bf c}^4 \rangle$, while the bulk permittivity $\boldsymbol\epsilon$ depends exclusively on $\langle {\bf c}^2 \rangle$.
To overcome this, two approaches can be taken.

### ➡️ Method 1 

Use the [IBOF closure](https://doi.org/10.1016/j.jnnfm.2005.11.005) of the Elmer/Ice flow model:
 
```python
import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) # L=4 is sufficient

### Determine <c^2> from radar-derived Delta lambda
l1 = 0   # lambda_1 = 0 (Gerber's approximation)
dl = 0.5 # Delta lambda = lambda_2 - lambda_1
a2 = np.diag([l1, l1+dl, 1-dl-2*l1]) # structure tensor <c^2> in eigenframe
a4 = sf.a4_IBOF(a2) # Elmer/Ice IBOF closure for <c^4> given <c^2>
```

### ➡️ Method 2

Use an empirical correlation for determining $\langle {\bf c}^4 \rangle$ given $\langle {\bf c}^2 \rangle$ if the CPO is approximately rotationally symmetric.

If the CPO symmetry axis is rotated into the vertical direction, $\langle {\bf c}^2 \rangle$ depends only on the normalized [spectral component](cpo-representation.md) $\hat{n}_2^0 = n_2^0/n_0^0:$

$$
\langle {\bf c}^2 \rangle = \frac{\bf I}{3} +  \frac{2\sqrt{5}}{15} \hat{n}_2^0
\left[\begin{matrix}
-1/2 & 0 & 0 \\ 
0  & -1/2  & 0 \\ 
0 & 0 & 1
\end{matrix}\right]
,
$$

and $\langle {\bf c}^4 \rangle$ only on $\hat{n}_2^0$ and $\hat{n}_4^0 = n_4^0/n_0^0$ (not shown).
The figure below shows the empirical correlation between these two components based on ice-core samples.

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/research/state-space/ice/state-space-empcorr.png){: style="width:540px"}

Thus, if $\hat{n}_2^0$ is extracted from $\langle {\bf c}^2 \rangle$ in this frame, $\hat{n}_4^0$ can be derived and hence $\langle {\bf c}^4 \rangle$ reconstructed.
To pose the CPO in the original, unrotated eigenframe (${\bf m}_1, {\bf m}_2, {\bf z}$), the resulting expansion series is finally rotated back, allowing eigenenhancements to easily be calculated using *specfab*.

The following code demonstrates how to take each step with *specfab*:

```python
--8<-- "docs/snippets/gallery-radar-method2.py"
```

For reference, the following plots show the different CPOs at each step for $\Delta\lambda=0.5$ and $\Delta\lambda=1$, respectively:

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/radar-PP-figs/code-example-output-dl0.5.png){: style="width:600px"}

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/radar-PP-figs/code-example-output-dl1.0.png){: style="width:600px"}


