## specfab

A spectral orientation fabric model for polycrystalline materials that:
- Can model lattice (c-axis) rotation, discontinuous dynamic recrystallization, and rotation recrystallization.
- Can calculate fabric-induced directional viscosities (enhancement factors) using Sachs/Taylor homogenizations.
- Can calculate elastic P- and S-wave velocities using Voigt/Reuss homogenizations.
- Contains expressions for both forward and inverse orthotropic and tranversely isotropic rheologies.

![image](demo/cube-crush-animation/cube-crush.gif)

## Examples
See the [specfab Wiki](https://github.com/nicholasmr/specfab/wiki)

## Install

| Enviroment | How to |
| :--- | :--- |
| Python | A pre-compiled module exists for Linux:<br>`pip3 install numpy --upgrade && pip3 install specfabpy` |
| Compile Python module |- For a local-only install, run `make specfabpy` in `/src` (requires LAPACK and BLAS) <br>- To install for general use (in other folders), run `make python` in `/src`. Note that if you do not have write permissions for your python installation, you can instead run `make specfabpy; python setup.py install --user`|
| Fortran | The Fortran module is built by running `make specfab.o` |
| Elmer/Ice Interface | To interface with Elmer/Ice, you need a shared version of the libraries (built with the same Fortran compiler as Elmer). If needed, change the compiler in `src/Makefile`, then run `make libspecfab.so` |

## Documentation

| Component | Reference |
| :--- | :--- |
| Lattice rotation | [Rathmann et al. (2021)](https://doi.org/10.1017/jog.2020.117) |
| Discontinuous dynamic recrystallization | [Rathmann and Lilien (2021)](https://doi.org/10.1017/jog.2021.88) |
| Orthotropic bulk rheologies | [Rathmann and Lilien (2022)](https://doi.org/10.1017/jog.2022.33) |

See also the [Wiki](https://github.com/nicholasmr/specfab/wiki)

## Q&A

| **Q** | **A** |
| :--- | :--- |
| What $L$ are possible? | Any $4\leq L\leq 20$. If higher $L$ are required: <br>1. `cd src/include && python3 make_gaunt_coefs.py L` (replacing `L`) <br>2. `make clean && make` |
