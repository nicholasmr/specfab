# Idealized CPOs

If concerned with the distribution of a *single* crystallographic axis, three types of idealized CPO states can be said to exist:

* **Unidirectional CPO**: all axes are perfectly aligned, i.e. perfect single maximum.
* **Planar CPO**: all axes are perfectly distributed on a plane, i.e. a great circle on $S^2$.
* **Circle CPO**: all axes are perfectly distributed on a small circle on $S^2$.

Each of these can be expanded as a spherical harmonics series by using the sifting property of the delta function $\delta({\hat {\bf r}})$.

## Unidirectional

Consider the case where slip-plane normals are perfectly aligned with ${\hat {\bf r}}_0$ such that $n({\hat {\bf r}}) = \delta({\hat {\bf r}}- {\hat {\bf r}}_0)$ where ${\hat {\bf r}}(\theta,\phi)$ is an arbitrary radial unit vector.
The corresponding expansion coefficients follow from the evaluating overlap integral 

$$
n_l^m 
= \int_{S^2} \delta(\hat{{\bf r}}-{\hat {\bf r}}_0) (Y_l^m(\hat{{\bf r}}))^* \,\mathrm{d}\Omega
= (Y_l^m({\hat {\bf r}}_0))^*
.
$$

In the figure below, the resulting unidirectional distribution is shown (rightmost inset), where the white area represents the subspace of admissible CPOs when expressed in terms of the normalized coefficients of lowest order: $\hat{n}_2^0 = n_2^0/n_0^0$ and $\hat{n}_4^0 = n_4^0/n_0^0$.

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/research/state-space/ice/state-space-ideal.png#center){: style="width:570px"}

The code below demonstrates how to generate the distribution with specfab.

```python
--8<-- "docs/snippets/idealized-cpo.py"
```

## Planar and circle

Planar and circle distributions follow from averaging the delta function over a desired co-latitude $\theta$ &mdash; i.e., the co-latitude where $n(\hat{{\bf r}})$ should be sharply defined &mdash; in which case all zonal structure vanishes ($m\neq 0$ components vanish) and we are left with

$$
n_l^m(\theta) = 
\begin{cases}
Y_l^0(\theta, \phi=0) \qquad\text{if}\quad m=0\\
0 \qquad\qquad\qquad\quad \text{if} \quad m\neq 0
\end{cases}
.
$$

Here, ${\hat {\bf r}}_0$ is to be understood as the rotational symmetry axis of $n(\hat{{\bf r}})$, and the co-latitude is defined w.r.t. ${\hat {\bf r}}_0$, not $\hat{{\bf z}}$.

The above figure also shows the resulting $n(\hat{{\bf r}})$ for different $\theta$, calculated using the same code as above but with a nonzero `colat`.

