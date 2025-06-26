# Ice viscous anisotropy

/// html | div[style='float: left; width: 30%; text-align: center;']
![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/monoice-viscous.png){: style="width:150px"}
///
/// html | div[style='float: right; width: 70%;']
If ice grains are treated as transversely isotropic, the rheology of a single grain can be modeled as a [transversely isotropic power law](constitutive-viscoplastic.md).
This requires specifying the grain eigenenhancements $E_{cc}'$ and $E_{ca}'$ and the power law exponent $n'$. 
///
/// html | div[style='clear: both;']
///

The grain parameters proposed by [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88) assume linear-viscous behavior of single crystals ($n'=1$) and promote the activation of basal glide by making that slip system soft compared to other systems: $E_{ca}' > 1$, whereas $E_{cc}'=1$. 
This reduces the problem to picking $E_{ca}'$ and $\alpha$ (Taylor&mdash;Sachs homogenization weight), which [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88) determined by requiring that deformation tests on strong single-maximum CPOs (aligned grains) are approximately reproduced; that is, $E_{mt}=10$ and $E_{mm}=0.01$.

The effect of choosing different $E_{ca}'$ and $\alpha$ (left panel) on the eigenenhancements of different CPO states (right panel) is shown below for combinations of $E_{ca}'$ and $\alpha$ that fulfill $E_{mt}=10$ given a unidirectional CPO. 

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/research/calibrate-Eij/ice/Eij-state-space/Eij-state-space.gif){: style="width:750px"}

Clearly, there is a trade-off between how shear enhanced ($E_{mt}$) and how hard for axial compression ($E_{mm}$) the model allows a unidirectional CPO to be.
    
### 📝 Code example

The below code shows how to calculate $E_{ij}$ given ${\bf a}^{(4)}$ (or the state vector $\bf s$) and a set of grain parameters.

```python
--8<-- "docs/snippets/Eij-LTS-ice.py"
```

