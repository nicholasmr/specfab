# Elastic wave propagation

Plane waves are a well-known solution to the Cauchy-Navier equation

$$
\nabla\cdot {\boldsymbol \sigma} = \rho \frac{\partial^2 {\bf d}}{\partial t^2},
$$

where ${\boldsymbol \sigma}$ is the bulk stress tensor, ${\bf d}$ is the displacement field, and $\rho$ the mass density.
Substituting ${\bf d}$ for the plane wave solution, ${\bf d} = {\bf d}_0 \exp[i({\bf k}\cdot {\bf x} - \omega t)]$, the problem reduces to

$$
(k^2\hat{{\bf Q}} - \omega^2 \rho {\bf I}) \cdot {\bf d} = {\bf 0},
$$

where $\hat{{\bf Q}}(\hat{{\bf k}})$ is the normalized [Christoffel (acoustic) tensor](https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Elasticity/Elasticity.htm), the form of which depends on the bulk [constitutive equation](constitutive-elastic.md) substituted for ${\boldsymbol \sigma}({\bf d})$.
The above equation implies that $\det( k^2\hat{{\bf Q}} - \omega^2 \rho {\bf I} ) = 0$. 
Hence, the eigenvalues and eigenvectors of $\hat{{\bf Q}}/\rho$ are the permitted wave velocities squared $V^2 = {\omega^2}/{k^2}$ and wave polarization, respectively.

## Homogenization

The problem may be closed by approximating $\hat{{\bf Q}}$ as the grain-averaged acoustic tensor, subject to a  linear combination of the Voigt and Reuss homogenization schemes:

$$
\hat{{\bf Q}} = (1-\alpha) \langle\hat{{\bf Q}}'_{\mathrm{Reuss}}\rangle + \alpha \langle \hat{{\bf Q}}'_{\mathrm{Voigt}}\rangle
.
$$    

The Voigt scheme ($\alpha=1$) assumes the strain field is homogeneous over the polycrystal scale, whereas the Reuss scheme ($\alpha=0$) assumes the stress is homogeneous.

In these homogenizations, grains are therefore assumed interactionless and the bulk elastic behaviour is therefore simply the grain-orientation-averaged elastic behaviour subject to either homogeneous stress or strain assumptions over the polycrystal scale.

<!--    
Elastic P and S plan-wave velocities, permitted in a polycrystal with an arbitrary CPO, can be determined for any linear combination of the Voigt and Reuss homogenization schemes by solving for the eigenvalues of the [acoustic tensor](https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Elasticity/Elasticity.htm):
-->
    
!!! warning "Grain parameters" 
    The grain elastic parameters used for homogenization should be understood as the *effective* polycrystal values needed to reproduce experimental results; they are not the values derived from experiments on single crystals.

- - -

## Glacier ice

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/monoice-elastic.png){: style="width:150px"}

If ice grains are approximated as transversely isotropic, their elastic behaviour can be modelled using the [transversely isotropic elastic constitutive equation](constitutive-elastic.md).
This requires specifying the grain elastic parameters $\lambda'$, $\mu'$, $\hat{\lambda}'$, $\hat{\mu}'$, $\hat{\gamma}'$, and the Voigt&mdash;Reuss weight $\alpha$.

### 📝 Code example

```python
--8<-- "docs/snippets/waveprop-elastic-ice.py"
```

The below animation shows directional P- and S-wave velocities for a CPO evolving under [uniaxial compression](deformation-kinematics.md) along ${\hat {\bf z}}$, relative to an isotropic CPO, when subject to [lattice rotation](fabdyn-LROT.md).

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/cube-crush-animation/S2-maps/S2-vi.gif){: style="width:660px"}

- - -

## Olivine

 ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/monooli-elastic-mi.png){: style="width:200px"} 
 
If olivine grains are approximated as orthotropic, their elastic behaviour can be modelled using the [orthotropic elastic constitutive equation](constitutive-elastic.md).
This requires specifying the grain elastic parameters $\lambda_{ij}'$, $\mu_{i}'$, and the Voigt&mdash;Reuss weight $\alpha$.

### 📝 Code example

```python
--8<-- "docs/snippets/waveprop-elastic-olivine.py"
```

