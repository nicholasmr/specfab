# Viscoplastic constitutive equations

Anisotropic power-law rheologies are supported in both forward and inverse (or reverse) form.

<!--
* Forward form: ${\bf D}({\bf S})$ 
* Inverse form: ${\bf S}({\bf D})$ 

where $\bf{D}$ and $\bf{S}$ are the strain-rate and deviatoric stress tensors, respectively.
-->

!!! note "Eigenenhancements"
    Anisotropic viscosities/fluidities are prescribed in terms of logitudinal and shear strain-rate enhancement factors w.r.t material symmetry axes, termed [eigenenhancements](enhancements-strainrate.md).

## Transversely isotropic

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous-bulk.png){: style="width:140px"} 

### `D = sf.rheo_fwd_tranisotropic(S, A, n, m, Emm, Emt)`

### `S = sf.rheo_rev_tranisotropic(D, A, n, m, Emm, Emt)`

| Arguments | Type |
| --- | --- |
| `S`, `D` | Deviatoric-stress and strain-rate tensor (3x3) |
| `A`, `n` | Flow-rate factor $A$ and power-law exponent $n$  |
| `m` | Rotational symmetry axis $\bf{m}$  |
| `Eij` | Tuple of eigenenhancements `(Emm, Emt)`|


## Orthotropic

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous-bulk.png){: style="width:250px"} 

### `D = sf.rheo_fwd_orthotropic(S, A, n, m1,m2,m3, Eij)`

### `S = sf.rheo_rev_orthotropic(D, A, n, m1,m2,m3, Eij)`

| Arguments | Type |
| --- | --- |
| `S`, `D` | Deviatoric-stress and strain-rate tensor (3x3) |
| `A`, `n` | Flow-rate factor $A$ and power-law exponent $n$  |
| `m1`, `m2`, `m3` | Reflection symmetry axes $\bf{m}_1$, $\bf{m}_2$, and $\bf{m}_3$  |
| `Eij` | Tuple of eigenenhancements `(E11,E22,E33,E23,E13,E12)` |

