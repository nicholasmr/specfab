# Continous dynamic recrystallization (CDRX) 

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/iceproc-CDRX.png){: style="width:120px"}

Polygonization (rotation recrystallization, CDRX) accounts for the division of grains along internal sub-grain boundaries resulting from local strain incompatibilities. 
In effect, CDRX reduces the average grain size upon grain division but does not necessarily change the CPO much ([Alley, 1992](https://doi.org/10.3189/S0022143000003658)). 

Following [Gödert (2003)](https://doi.org/10.1007/s001610050095), CDRX can be modeled by approximating this effect as a Laplacian diffusive process on $S^2$:

$$
\frac{\mathrm{D} n}{\mathrm{D} t} = \Lambda\nabla^2 n  ,
$$

where $\Lambda$ is the CDRX rate-factor magnitude that depends on temperature, stress, strain-rate, etc. ([Richards et al., 2021](https://doi.org/10.1016/j.epsl.2020.116718)).

!!! warning "Limited use"

    Notice that this model is currently only relevant to slip-system normals. 
    The model is therefore not yet useful for e.g. olivine. 

### Matrix model 

The corresponding effect on the continuous distribution function is 

$$
\frac{\mathrm{D} n}{\mathrm{D} t} = \Lambda\nabla^2 n 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}}{\mathrm{D} t} = {\bf M_{\mathrm{CDRX}}} \cdot {\bf s} .
$$

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/fabric-evolution/animation-CDRX.gif){: style="width:610px"}

### Code example 

```python
--8<-- "docs/snippets/fabdyn-CDRX.py"
```

