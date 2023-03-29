# Elastic constitutive equations

Linear elastic constituve equations are supported in both forward and inverse (or reverse) form:

- Forward form ${\bf E}({\bf S})$ 
- Inverse form ${\bf S}({\bf E})$ 

where $\bf{E}$ and $\bf{S}$ are the strain and stress tensors, respectively.

## Transversely isotropic

| Symmetries |  Stiffness matrix ${\bf C}$ for $\bf{m}=\hat{\bf{z}}$ | 
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/tranisotropic/tranisotropic-elastic-bulk.png){: style="width:120px"} | $\begin{bmatrix}\gamma & \lambda & \hat{\lambda}\lambda & 0 & 0 & 0 \\\lambda & \gamma & \hat{\lambda}\lambda & 0 & 0 & 0 \\\hat{\lambda}\lambda & \hat{\lambda}\lambda & \hat{\gamma}\gamma & 0 & 0 & 0 \\0&0&0& \hat{\mu}\mu & 0 & 0\\0&0&0& 0 & \hat{\mu}\mu & 0\\0&0&0& 0 & 0 & \mu\\\end{bmatrix}$       | 

### `E = sf.elas_fwd_tranisotropic(S, lam, mu, Elam, Emu, Egam, m)`

### `S = sf.elas_rev_tranisotropic(E, lam, mu, Elam, Emu, Egam, m)`

where 

| Arguments | Type |
| --- | --- |
| `S`, `E` | stress and strain tensor (3x3) |
| `lam`, `mu` | Isotropic Lamé parameters $\lambda$ and $\mu$  |
| `Elam`, `Emu`, `Egam` | Anisotropic enhancement factors $\hat{\lambda}$, $\hat{\mu}$, and $\hat{\gamma}$  |

!!! note
    P-wave modulus is $\gamma \equiv \lambda + 2\mu$ (i.e. not a free parameter).
    
!!! tip "Tip: convert from $C_{ij}$ to Lamé parameters"

    ```
    lam,mu,Elam,Emu,Egam = sf.Cij_to_Lame_tranisotropic(C11,C33,C55,C12,C13) 
    ```

## Orthotropic

| Symmetries |  Stiffness matrix ${\bf C}$ for $({\bf m}_1,{\bf m}_2,{\bf m}_3)=(\hat{{\bf x}},\hat{{\bf y}},\hat{{\bf z}})$ | 
| :-: | :-: |
| ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/orthotropic/orthotropic-elastic-bulk.png){: style="width:250px"} | $\small\begin{bmatrix} \lambda_{11} + 2\mu_1 & \lambda_{12} & \lambda_{13} & 0 & 0 & 0 \\ \lambda_{12} & \lambda_{22} + 2\mu_2 & \lambda_{23} & 0 & 0 & 0 \\ \lambda_{13 }& \lambda_{23} & \lambda_{33} + 2\mu_3 & 0 & 0 & 0 \\ 0 & 0 & 0 & \dfrac{\mu_2+\mu_3}{2} & 0 & 0 \\ 0 & 0 & 0 & 0 & \dfrac{\mu_3+\mu_1}{2} & 0 \\ 0 & 0 & 0 & 0 & 0 & \dfrac{\mu_1+\mu_2}{2} \end{bmatrix}$ | 

Not yet available

<!--

### `E = sf.elas_fwd_orthotropic(S, lam, mu, Elam, Emu, Egam, m)`

### `S = sf.elas_rev_orthotropic(E, lam, mu, Elam, Emu, Egam, m)`

where 

| Arguments | Type |
| --- | --- |
| `S`, `E` | stress and strain tensor (3x3) |
| `lam`, `mu` | Isotropic Lamé parameters $\lambda$ and $\mu$  |
| `Elam`, `Emu`, `Egam` | Anisotropic enhancement factors $\hat{\lambda}$, $\hat{\mu}$, and $\hat{\gamma}$  |

-->
