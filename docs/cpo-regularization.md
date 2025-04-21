# Regularization

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/demo/cube-crush-animation/regularization/regularization.gif){: style="width:700px"}

As $n(\theta,\phi)$ and $b(\theta,\phi)$ becomes anisotropic due to CPO processes, the coefficients $n_l^m$ and $b_l^m$ associated with high wavenumber modes (large $l$ and $m$, and thus small-scale structure) must increase in magnitude relative to the low wavenumber coefficients (small $l$ and $m$). 

One way to visualize this is by the [angular power spectrum](https://en.wikipedia.org/wiki/Spherical_harmonics#Spectrum_analysis) (and similarly for $b$)

$$ 
S(l) = \frac{1}{2l + 1} \sum_{m=-l}^l \left\vert n_l^m \right\vert^2 ,
$$

which grows with time. 
In the animation above, the left-hand panel shows how the power spectrum evolves under lattice rotation (unconfined vertical compression) compared to the end-member case of a delta function (dashed line).

If the expansion series is truncated at $l=L$, then $l{\gt}L$ modes cannot evolve, and the truncated solution will reach an unphysical quasi-steady state. To prevent this, regularization must be introduced.
Specfab uses Laplacian hyper diffusion ($k>1$) as regularization in $S^2$

$$ 
\frac{\mathrm{D} n_l^m}{\mathrm{D} t} = \sqrt{\frac{{\bf D}:{\bf D}}{2}} {\nu}[l(l+1)]^{k} n_l^m 
\quad\Longrightarrow\quad
\frac{\mathrm{D} {\bf s}}{\mathrm{D} t} = {\bf M_{\mathrm{REG}}} \cdot {\bf s} 
,
$$

that can be added to the fabric evolution operator ${\bf M}$ as follows:
```python
M += sf.M_REG(nlm, D)
```
This allows the growth of high wavenumber modes to be disproportionately damped (green line compared to red line in animation above).
The dependency on the strain-rate tensor guarantees the same behavior irrespective of strain rate magnitude. 


!!! note "Rule-of-thumb"
    As a rule-of-thumb, regularization affects the highest and next-highest modes $l{\geq}L-2$ and can therefore not be expected to evolve freely. 
    This, in turn, means that [structure tensors](cpo-representation.md) `a2` and `a4`, and hence [calculated enhancement factors](enhancements-strainrate.md), might be affected by regularization unless $L{\geq}8$ is chosen. 

