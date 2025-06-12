# Lattice rotation

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/slip-plane/plastic-spin.png){: style="width:350px"}

When subject to simple shear, the orientation of easy slip systems in polycrystalline materials like ice or olivine tend towards aligning with the bulk shear-plane system. 
That is, the $n$-axis ($c$-axis in ice) tends towards to bulk shear plane normal, and the $b$-axis ($a$-axis in ice) towards to bulk shear direction.
Thus, if grains are perfectly aligned with the bulk shear system, their orientation should be unaffected by any further shearing, but be in steady state.
Clearly, slip systems do therefore not simply co-rotate with the bulk continuum spin (${\boldsymbol \omega}$) like passive material line elements embedded in a flow field, i.e. 
${\bf \dot{n}} \neq {\boldsymbol \omega} \cdot {\bf n}$.
Rather, slip systems must be subject to an additional contribution &mdash; a plastic spin ${\boldsymbol\omega}'$ &mdash; such that the bulk spin is exactly counteracted to achieve steady state if favorably aligned:

<!-- \quad\text{for ${\bf b}$&ndash;${\bf n}$ shear}. -->
$$
{\bf \dot{n}} = ({\boldsymbol\omega} + {\boldsymbol\omega}'_{n}) \cdot {\bf n} = {\bf 0} , 
$$

$$
{\bf \dot{b}} = ({\boldsymbol\omega} + {\boldsymbol\omega}'_{b}) \cdot {\bf b} = {\bf 0} .
$$

More precisely, the crystallographic axes reorient themselves in response to both the bulk continuum spin and a plastic spin that is supposed to represent the crystallographic spin needed to accommodate strain compatibility between grains that otherwise preferentially deform by easy slip (i.e., basal slip for ice).

## Directors method

[Aravas and Aifantis (1991)](https://doi.org/10.1016/0749-6419(91)90028-W) and [Aravas (1994)](https://www.doi.org/10.1088/0965-0393/2/3A/005) (among others) proposed a particularly simple model for the functional form of ${\boldsymbol \omega}'$, the so-called directors method.

For a constant rate of shear deformation ($1/\tau$) aligned with the ${\bf b}$&mdash;${\bf n}$ system, the deformation gradient is ${\bf F} = {\bf I} + (t/\tau) {\bf b}{\bf n}$ and therefore 

$$
\nabla {\bf u} = \frac{1}{\tau} {\bf b}{\bf n}
\quad \Longrightarrow\quad
{\dot{\boldsymbol\epsilon}} = \frac{1}{2\tau} ({\bf b}{\bf n} + {\bf n}{\bf b}),
\quad
{\boldsymbol \omega} = \frac{1}{2\tau} ({\bf b}{\bf n} - {\bf n}{\bf b}) .
$$

Since ${\boldsymbol \omega}' = -{\boldsymbol \omega}$ is required in steady state, it follows from eliminating $1/(2\tau)$ by calculating 
${\dot{\boldsymbol\epsilon}} \cdot {\bf n}^2$, 
${\bf n}^2 \cdot {\dot{\boldsymbol\epsilon}}$, 
${\dot{\boldsymbol\epsilon}} \cdot {\bf b}^2$, and 
${\bf b}^2 \cdot {\dot{\boldsymbol\epsilon}}$, 
that

$$
{\boldsymbol \omega}'_{n} = +{\bf n}^2 \cdot {\dot{\boldsymbol\epsilon}} - {\dot{\boldsymbol\epsilon}} \cdot {\bf n}^2 ,
\\
{\boldsymbol \omega}'_{b} = -{\bf b}^2 \cdot {\dot{\boldsymbol\epsilon}} + {\dot{\boldsymbol\epsilon}} \cdot {\bf b}^2 .
$$

Indeed, this result agrees with representation theorems for isotropic functions (Wang, 1969), stating that an antisymmetric tensor-valued function of a symmetric tensor (${\dot{\boldsymbol\epsilon}}$) and a vector (${\hat {\bf r}}$) is to lowest order given by

$$
{\boldsymbol \omega}' = 
\iota({\hat {\bf r}}^2\cdot{\dot{\boldsymbol\epsilon}} - {\dot{\boldsymbol\epsilon}}\cdot{\hat {\bf r}}^2)
.
$$

To be consistent with the above, $\iota = +1$ for ${\hat {\bf r}}={\bf n}$ and $\iota = -1$ for ${\hat {\bf r}}={\bf b}$.

${\boldsymbol \omega}'_{n}$ and ${\boldsymbol \omega}'_{b}$ are then generally taken to suffice for other deformation kinematics, too.

!!! note "Glacier ice"

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/plastic-spin.png){: style="width:350px"}

    The predicted normalized $c$-axis velocity field for glacier ice where ${\bf n} = {\bf c}$ is show below for three different deformation kinematics:

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/fabric-dynamics/fabdyn-LROT.png#center){: style="width:470px"}

### Matrix model 

The corresponding effect on the continuous distribution functions is modelled as a conservative advection process in orientation space $S^2$:

$$ 
\frac{\mathrm{D} n}{\mathrm{D} t} = -\nabla_{S^2}\cdot(n{\bf \dot{n}}) 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}_n}{\mathrm{D} t} = {\bf M_{\mathrm{LROT}}}(\iota=+1) \cdot {\bf s}_n,
$$

$$ 
\frac{\mathrm{D} b}{\mathrm{D} t} = -\nabla_{S^2}\cdot(b{\bf \dot{b}}) 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}_b}{\mathrm{D} t} = {\bf M_{\mathrm{LROT}}}(\iota=-1) \cdot {\bf s}_b,
$$

where ${\bf M_{\mathrm{LROT}}}(\iota)$ is given analytically in [Rathmann et al. (2021)](https://doi.org/10.1017/jog.2020.117).

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/fabric-evolution/animation-LROT.gif){: style="width:610px"}


### Code example

```python
--8<-- "docs/snippets/fabdyn-LROT.py"
```

