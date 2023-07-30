# Idealized CPOs

Three types of idealized CPO states can be said to exist:

* **Unidirectional CPO**: crystallographic axes are perfectly aligned, i.e. perfect single maximum.
* **Planar CPO**: crystallographic axes are perfectly distributed on a plane, i.e. a great circle on $S^2$.
* **Circle CPO**: crystallographic axes are perfectly distributed on a small circle on $S^2$.

Each of these can be expanded in terms of spherical harmonics by using the sifting property of the delta function, $\delta({\hat {\bf r}})$.

## Unidirectional distribution

Consider the case where grains are perfectly aligned with ${{\bf m}}$ such that $n({\hat {\bf r}}) = \delta({\hat {\bf r}}-{{\bf m}})$.
The corresponding expansion coefficients follow from the usual overlap integral:

$$
n_l^m 
= \int_{S^2} \delta(\hat{{\bf r}}-{{\bf m}}) (Y_l^m(\hat{{\bf r}}))^* \,\mathrm{d}\Omega
= (Y_l^m({{\bf m}}))^*
.
$$

In the figure below, the resulting unidirectional distribution is shown (right-most inset), where the white area represents the lowest-order CPO state space given by the normalized coefficients $\hat{n}_2^0 = n_2^0/n_0^0$ and $\hat{n}_4^0 = n_4^0/n_0^0$:

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/state-space-validation/ice/state-space-ideal.png#center){: style="width:570px"}

The code below demonstrates how to generate the distribution with specfab.

```python
import numpy as np
from specfabpy import specfabpy as sf
L = 8
lm, nlm_len = sf.init(L) 

m = [0,0,1] # symmetry axis of distribution 
colat = 0 # 0 = unidirectional distribution, pi/2 = planar distribution, and anything in between is a small circle distribution
nlm = sf.nlm_ideal(m, colat, L) # note: only l<=12 coefs are determined even if L>12
```

## Planar and circle distribution

Planar and circle distributions follow from averaging the delta function over a desired co-latitude $\theta$ &mdash; i.e. the co-latitude where $n(\hat{{\bf r}})$ should be sharply defined &mdash; in which case all zonal structure vanishes ($m\neq 0$ components vanish) and we are left with

$$
n_l^m(\theta) = 
\begin{cases}
Y_l^0(\theta, \phi=0) \qquad\text{if}\quad m=0\\
0 \qquad\qquad\qquad\quad \text{if} \quad m\neq 0
\end{cases}
.
$$

Here, ${{\bf m}}$ is to be understood as the rotational symmetry axis of $n(\hat{{\bf r}})$, and the co-latitude is defined w.r.t. ${{\bf m}}$ (not $\hat{{\bf z}}$).

The above figure also shows the resulting $n(\hat{{\bf r}})$ for different $\theta$, calculated using the same code as above but for nonzero `colat`.

