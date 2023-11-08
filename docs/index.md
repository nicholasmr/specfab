# specfab documentation

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/logo-square.jpg){: style="width:300px"}

Spectral CPO model of polycrystalline materials that can:

- Model lattice rotation, discontinuous DRX, and rotation/continuous DRX.
- Calculate CPO-induced viscous anisotropies using Sachs/Taylor homogenizations.
- Calculate elastic P- and S-wave velocities using Voigt/Reuss homogenizations.
- Provide expressions for forward+inverse orthotropic and transversely isotropic rheologies.
- Convert between structure tensors and spectral expansions coefficients.
- Be integrated with finite-element codes such as Elmer and FEniCS.

By Nicholas M. Rathmann and David A. Lilien

## Glacier ice demo

![](https://github.com/nicholasmr/specfab/raw/main/demo/cube-crush-animation/Eij-trajectory/Eij-trajectory.gif){: style="width:700px"}

<!-- ![](https://github.com/nicholasmr/specfab/raw/main/images/tranisotropic/parcel-animation/tranisotropic-parcel-animation.gif){: style="width:550px"} -->

## Install

Source code is [available here](https://github.com/nicholasmr/specfab).

| Environment | How to install |
| :--- | :--- |
| Python in Linux | - PyPI package: `pip3 install numpy --upgrade && pip3 install specfabpy`<br> - Compile yourself: `cd src && make python` |
| Python in Windows/Mac | You will have to compile specfab yourself. |
| Fortran | Run `cd src && make specfab.o` |
| Elmer/Ice Interface | Compile shared library by running `cd src && make libspecfab.so` |

Libraries required: BLAS, LAPACK

## Initialize 

Initialize `specfab` by running

```python
import numpy as np
from specfabpy import specfab as sf

lm, nlm_len = sf.init(10) # L=10 truncation is sufficient for many cases

nlm = np.zeros(nlm_len, dtype=np.complex64) # vector of harmonic expansion coefficients
nlm[0] = 1/np.sqrt(4*np.pi) # normalized isotropic distribution
```

where

| Variable | Interpretation |
| --- | --- |
| `nlm_len` | Number of expansion coefficients for expansion series truncated at $l=L$ |
| `nlm`     | Vector of complex-valued expansion coefficients (state vector) |
| `lm`      | Vector of degree and order integers (`l`,`m`) associated with each entry in `nlm` |

