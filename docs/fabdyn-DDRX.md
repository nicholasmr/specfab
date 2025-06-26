## Discontinous dynamic recrystallization (DDRX) 

/// html | div[style='float: left; width: 25%; text-align: center;']
![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/iceproc-DDRX.png){: style="width:120px"}
///
/// html | div[style='float: right; width: 75%;']
Following [Placidi and others (2010)](https://doi.org/10.1007/s00161-009-0126-0), DDRX is modeled as a spontaneous mass decay&mdash;production process in orientation space $S^2$, intended to represent the combined effect of nucleation and grain boundary migration. 
That is, mass is spontaneously exchanged between grains with different orientations depending on the local stress state, strain rate, and temperature, in a statistical sense. 
///
/// html | div[style='clear: both;']
///


Following [Placidi and others (2010)](https://doi.org/10.1007/s00161-009-0126-0), DDRX is modeled as a spontaneous mass decay&mdash;production process in orientation space $S^2$, intended to represent the combined effect of nucleation and grain boundary migration. 
That is, mass is spontaneously exchanged between grains with different orientations depending on the local stress state, strain rate, and temperature, in a statistical sense. 

The decay&mdash;production rate is defined as

$$
\Gamma = \Gamma_0\left(D- {\langle} D {\rangle}\right) 
$$

where the  prefactor $\Gamma_0$ accounts for the preferential (Arrhenius) activation at warm temperatures and the effect of strain-rate magnitude, defined as
\begin{align}
\Gamma_0 = \sqrt{\frac{{\dot{\boldsymbol\epsilon}}:{\dot{\boldsymbol\epsilon}}}{2}} A_{\Gamma}\exp(-Q_{\Gamma}/RT)
.
\end{align}
Here, ${\dot{\boldsymbol\epsilon}}$ is the strain-rate tensor, $R$ is the gas constant, $T$ is the temperature, and $A_{\Gamma}$ and $Q_{\Gamma}$ have been calibration by [Richards et al. (2021)](https://doi.org/10.1016/j.epsl.2020.116718) and [Lilien et al. (2023)](https://doi.org/10.1017/jog.2023.78). 

The deformability $D$ is the square of the basal-plane resolved shear stress

$$
\tau^2_\mathrm{RSS} = ({\boldsymbol\tau}\cdot{\boldsymbol\tau}):{\bf n}^2 - {\boldsymbol\tau}:{\bf n}^4:{\boldsymbol\tau},
$$

relative to the average value of an isotropic fabric:

$$
D = \frac{\tau^2_\mathrm{RSS}}{\langle \tau^2_\mathrm{RSS} \rangle_\mathrm{iso}} = \frac{({\boldsymbol\tau}\cdot{\boldsymbol\tau}):{\bf n}^2 - {\boldsymbol\tau}:{\bf n}^4:{\boldsymbol\tau}}{{\boldsymbol\tau}:{\boldsymbol\tau}/5},
$$

where ${\bf n}$ is an arbitrary slip-plane normal ($c$-axis for ice). 
Because $D$ is largest for grains with an orientation favorable to easy glide, mass is spontaneously created/added to grains with such preferred orientations (in a statistical sense).
Conversely, mass spontaneously decays if $D<{\langle} D {\rangle}$, corresponding to grains with an unfavorable orientation being consumed by grains with a more favorable orientation to basal glide. 
Here, ${\langle} D {\rangle}$ is the grain-average deformability of the polycrystal, and the total mass of the polycrystal is conserved since $\Gamma$ is zero on average, $\langle\Gamma\rangle=\Gamma_0\langle D-\langle D\rangle\rangle=0$ (equal amounts of grain mass are created and destroyed per time).

The primary driver for grain boundary migration is generally regarded to be the contrast in dislocation density (stored strain energy) between grains. 
In this sense, $1/D$ can be understood as a parameterization of the dislocation density: the larger $1/D$ is, the more dislocations a grain has, and the more likely it is to be consumed (decay).

!!! warning "Nonlinear process"

    Since average deformability $\langle D\rangle$ depends on the instantaneous CPO state &mdash; specifically, the [structure tensors](cpo-representation.md) $\langle{\bf c}^2\rangle$ and $\langle{\bf c}^4\rangle$ &mdash; this crystal process is nonlinear (renders a nonlinear matrix problem below).

!!! warning "Limited use"

    Notice that this model is currently only relevant to slip-system normals. 
    The model is therefore not yet useful for e.g. olivine. 

!!! note "Glacier ice"

    The normalized decay&mdash;production rate is shown below for three different stress states in the case of glacier ice where ${\bf n} = {\bf c}$:

    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/fabric-dynamics/fabdyn-DDRX.png#center){: style="width:470px"}

### Matrix model 

The corresponding effect on the continuous distribution function is 

$$ 
\frac{\mathrm{D} n}{\mathrm{D} t} = \Gamma n 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}_n}{\mathrm{D} t} = {\bf M_{\mathrm{DDRX}}} \cdot {\bf s}_n ,
$$

where ${\bf M_{\mathrm{DDRX}}}$ is given analytically in [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88).

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/fabric-evolution/animation-DDRX.gif){: style="width:610px"}

### 📝 Code example

```python
--8<-- "docs/snippets/fabdyn-DDRX.py"
```

