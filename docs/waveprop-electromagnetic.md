# Electromagnetic wave propagation

Plane waves are a well-known solution to Maxwell's equations in a non-conducting, source-free, anisotropic linear dielectric medium

$$
\nabla \times \nabla \times {\bf E} = -\mu {\boldsymbol \epsilon} \frac{\partial^2 {\bf E}}{\partial t^2},
$$

<!--Note this can also be written as $\nabla^2 {\bf E} = \mu {\boldsymbol \epsilon} {\partial^2 {\bf E}}/{\partial t^2}$.-->
where ${\bf E}$ is the electric field, ${\boldsymbol \epsilon}$ is the bulk dielectric permittivity tensor, and $\mu$ the bulk isotropic permeability of the medium.
Substituting ${\bf E}$ for the plane wave solution, ${\bf E} = {\bf E}_0 \exp[i({\bf k}\cdot {\bf x} - \omega t)]$, the problem reduces to

$$
({\bf K} + \omega^2 \mu {\boldsymbol \epsilon}) \cdot {\bf E} = {\bf 0},
$$

where ${\bf K} = {\bf k}\times{\bf k} \times$ is the [matrix representation](https://en.wikipedia.org/wiki/Cross_product#Alternative_ways_to_compute) of the twice-applied cross product.
The above equation implies that $\det( {\bf K} + \omega^2 \mu {\boldsymbol \epsilon} ) = 0$. Hence, the eigenvalues and eigenvectors of the normalized Christoffel tensor ${\boldsymbol \epsilon}^{-1} {\bf K}/(\mu k^2)$ are the permitted wave velocities squared $V^2 = {\omega^2}/{k^2}$ and wave polarization, respectively.

## Homogenization 

If wave lengths are much longer than the average grain size, the problem can be closed by approximating the bulk polycrystalline permittivity tensor by the grain-averaged permittivity tensor

$${\boldsymbol \epsilon} = \langle {\boldsymbol \epsilon}' \rangle,$$
    
constructed by averaging over all grain orientations (over the CPO).

- - -

## Glacier ice

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/monoice-dielectric.png){: style="width:150px"}

If ice grains are approximated as transversely isotropic w.r.t. the $c$-axis, the dielectric permittivity tensor of a single crystal can be written as
$$
{\boldsymbol \epsilon}' =  (2\epsilon_{a}' + \epsilon_{c}') \frac{\bf I}{3}
+ (\epsilon_{c}'-\epsilon_{a}') \left( {\bf c}^2 - \frac{\bf I}{3} \right),
$$
where $\epsilon_{c}'$ and $\epsilon_{a}'$ are the permittivities parallel and perpendicular to the $c$ axis.
In this case, the grain-averaged permitivity is simply 

$$
\langle {\boldsymbol \epsilon}'\rangle = (2\epsilon_{a}' + \epsilon_{c}') \frac{\bf I}{3}
+ (\epsilon_{c}'-\epsilon_{a}') \left( \langle {\bf c}^2 \rangle - \frac{\bf I}{3} \right)
,
$$

where $\langle {\bf c}^2 \rangle$ is the [second-order structure tensor](cpo-structuretensors.md).

### Code example

Experimental, bug reports are welcome.


```python
--8<-- "docs/snippets/waveprop-em-ice.py"
```

The below animation shows directional S-wave velocities for a CPO evolving under [$xz$ confined compression](deformation-kinematics.md), relative to an isotropic CPO, when subject to [lattice rotation](fabdyn-LROT.md).

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/cube-crush-animation/S2-maps/S2-vi-EM.gif){: style="width:460px"}

- - -

## Olivine

🚧 *Not yet supported.*
