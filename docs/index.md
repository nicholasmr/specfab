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

![](https://github.com/nicholasmr/specfab/raw/main/demo/cube-crush-animation/cube-crush.gif){: style="width:550px"}

## Install

Source code [available here](https://github.com/nicholasmr/specfab)

| Environment | How to |
| :--- | :--- |
| Python | A pre-compiled module exists for Linux:<br>`pip3 install numpy --upgrade && pip3 install specfabpy` |
| Compile Python module |- For a local-only install, run `make specfabpy` in `/src` (requires LAPACK and BLAS) <br>- To install for general use (in other folders), run `make python` in `/src`. Note that if you do not have write permissions for your python installation, you can instead run `make specfabpy; python setup.py install --user`|
| Fortran | The Fortran module is built by running `make specfab.o` |
| Elmer/Ice Interface | To interface with Elmer/Ice, you need a shared version of the libraries (built with the same Fortran compiler as Elmer). If needed, change the compiler in `src/Makefile`, then run `make libspecfab.so` |

## Initialize 

Initialize `specfabpy` by running

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(10) # L=10 truncation is sufficient for many cases
nlm = np.zeros(nlm_len, dtype=np.complex64)
nlm[0] = 1/np.sqrt(4*np.pi) # Normalized, isotropic distribution
```

where

| Variable | Interpretation |
| --- | --- |
| `nlm_len` | Number of expansion coefficients for expansion series truncated at $l=L$ |
| `nlm`     | Vector of complex-valued expansion coefficients (state vector) |
| `lm`      | Vector of degree and order integers (`l`,`m`) associated with each entry in `nlm` |

## Literature 

| Component | Reference |
| :--- | :--- |
| Lattice rotation | [Rathmann et al. (2021)](https://doi.org/10.1017/jog.2020.117) |
| Discontinuous dynamic recrystallization | [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88) |
| Orthotropic bulk rheologies | [Rathmann and Lilien (2022)](https://doi.org/10.1017/jog.2022.33) |
| Elastic wave velocities | [Rathmann et al. (2022)](https://doi.org/10.1098/rspa.2022.0574) |

