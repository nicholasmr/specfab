# Viscoplastic constitutive equations

Anisotropic power-law rheologies are supported in both forward and inverse (or reverse) form.

<!--
* Forward form: ${\bf D}({\bf S})$ 
* Inverse form: ${\bf S}({\bf D})$ 

where $\bf{D}$ and $\bf{S}$ are the strain-rate and deviatoric stress tensors, respectively.
-->

!!! note "Eigenenhancements"
    Anisotropic viscosities/fluidities are prescribed in terms of logitudinal and shear strain-rate enhancement factors w.r.t rheological symmetry axes, termed [eigenenhancements](enhancements-strainrate.md) ($E_{ij}$).

## Transversely isotropic

<!--
![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous-bulk.png){: style="width:140px"} 

### Forward rheology
--> 

| Rheolgical symmetries | Forward rheology |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-viscous-bulk.png){: style="width:140px"}  | $$ {\bf D} = \eta^{-1} \Big( {\bf S} - \lambda_1 ({\bf S}:{\bf M}){\bf I} + \lambda_2 ({\bf S}:{\bf M}){\bf M} + \lambda_3 ({\bf S}\cdot{\bf M} + {\bf M}\cdot{\bf S}) \Big) $$ $$ \eta^{-1} = A\Big( {\bf S}:{\bf S} + \lambda_2 ({\bf S}:{\bf M})^2 + 2\lambda_2 ({\bf S}^2:{\bf M}) \Big)^{(n-1)/2} $$ $$ {\bf M}={\bf m}^2$$|

<!--
$$
{\bf D} = \eta^{-1} \Big(
{\bf S} 
- \lambda_1 ({\bf S}:{\bf m}^2){\bf I} 
+ \lambda_2 ({\bf S}:{\bf m}^2){\bf m}^2
+ \lambda_3 ({\bf S}\cdot{\bf m}^2 + {\bf m}^2\cdot{\bf S})
\Big)
\\
\eta^{-1} = A\Big( {\bf S}:{\bf S} + \lambda_2 ({\bf S}:{\bf m}^2)^2 + 2\lambda_2 ({\bf S}^2:{\bf m}^2) \Big)^{(n-1)/2}
$$
-->

where the material parameters $\lambda_i$ depend on the eigenenhancements:

$$
\lambda_1 = \frac{E_{mm}^{2/(n+1)}-1}{2} ,\quad
\lambda_2 = \frac{3(E_{mm}^{2/(n+1)}-1) - 4(E_{mt}^{2/(n+1)}-1)}{2} ,\quad
\lambda_3 = E_{mt}^{2/(n+1)}-1
$$

### API

### `D = sf.rheo_fwd_tranisotropic(S, A, n, m, Eij)`

### `S = sf.rheo_rev_tranisotropic(D, A, n, m, Eij)`

| Arguments | Type |
| --- | --- |
| `S`, `D` | Deviatoric-stress and strain-rate tensor (3x3) |
| `A`, `n` | Flow-rate factor $A$ and power-law exponent $n$  |
| `m` | Rotational symmetry axis $\bf{m}$  |
| `Eij` | Tuple of eigenenhancements `(Emm, Emt)`|


## Orthotropic

<!--
![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous-bulk.png){: style="width:250px"} 

### Forward rheology
-->

| Rheolgical symmetries | Forward rheology |
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-viscous-bulk.png){: style="width:250px"} | ${\bf D} = \eta^{-1} \displaystyle\sum_{i=1}^{3} \Big[ \lambda_i ({\bf S}:{\bf M}_i){\bf M}_{i}  + \lambda_{i+3} ({\bf S}:{\bf M}_{i+3}) {\bf M}_{i+3} \Big]$ </br></br> $\eta^{-1} = A\left( \displaystyle\sum_{i=1}^3 \Big[ \lambda_i ({\bf S}:{\bf M}_{i})^2 + \lambda_{i+3} ({\bf S}:{\bf M}_{i+3})^2 \Big] \right)^{(n-1)/2}$ </br></br>  ${\bf M}_{i} = \dfrac{{\bf m}_{j_i} {\bf m}_{j_i} - {\bf m}_{k_i} {\bf m}_{k_i}}{2} ,\quad {\bf M}_{i+3} = \dfrac{{\bf m}_{j_i} {\bf m}_{k_i} + {\bf m}_{k_i} {\bf m}_{j_i}}{2},$ </br> </br> $(j_1, j_2, j_3) = (2,3,1),\quad (k_1, k_2, k_3) = (3,1,2)$|

<!--
$$
{\bf D} = \eta^{-1} \sum_{i=1}^{3} \left(
\lambda_i I_{i} \frac{ {\bf m}_{j_i}{\bf m}_{j_i} - {\bf m}_{k_i}{\bf m}_{k_i} }{2}
+
\lambda_{i+3} I_{i+3} \frac{ {\bf m}_{j_i}{\bf m}_{k_i} + {\bf m}_{k_i}{\bf m}_{j_i} }{2}
\right)
\\
\eta^{-1} = A\left( \sum_{i=1}^3 \left( \lambda_i I_i^2 + \lambda_{i+3} I_{i+3}^2 \right) \right)^{(n-1)/2}
$$
-->

where the material parameters $\lambda_i$ depend on the eigenenhancements:

$$
\lambda_i     = \frac{4}{3} \left( E_{j_i j_i}^{2/(n+1)} + E_{k_i k_i}^{2/(n+1)} - E_{i i}^{2/(n+1)} \right),\quad
\lambda_{i+3} = 2 E_{j_i k_i}^{2/(n+1)}
$$

### API

### `D = sf.rheo_fwd_orthotropic(S, A, n, m1,m2,m3, Eij)`

### `S = sf.rheo_rev_orthotropic(D, A, n, m1,m2,m3, Eij)`

| Arguments | Type |
| --- | --- |
| `S`, `D` | Deviatoric-stress and strain-rate tensor (3x3) |
| `A`, `n` | Flow-rate factor $A$ and power-law exponent $n$  |
| `m1`, `m2`, `m3` | Reflection symmetry axes $\bf{m}_1$, $\bf{m}_2$, and $\bf{m}_3$  |
| `Eij` | Tuple of eigenenhancements `(E11,E22,E33,E23,E13,E12)` |

