# Spectral orientation fabric model for polycrystals

<img src="https://raw.githubusercontent.com/nicholasmr/specfab/main/images/logo.jpg" width="200">

A spectral fabric model for polycrystalline materials (e.g. ice, olivine, etc.) that 
- Can model fabric processes such as kinematic lattice (c-axis) rotation (Svendsen and
Hutter, 1996), discontinuous dynamic recrystallization (Placidi and others, 2010), and rotation recrystallization (Gödert, 2003).
- Can calculate fabric-induced directional viscosities (directional strain-rate enhancement factors) using Sachs/Taylor homogenization schemes.
- Contains analytical expressions for both forward and inverse anisotropic rheologies such as the transversely isotropic rheology and the orthotropic rheology.

## Tutorials, python examples, and documentation
See the [specfab Wiki](https://github.com/nicholasmr/specfab/wiki)

## Demo of lattice rotation and resulting strain-rate enhancement factors
![image](demo/cube-crush-animation/cube-crush.gif)

## Installation
### Python module `specfabpy`
- If you are running Linux, a pre-compiled version is available: `pip3 install numpy --upgrade && pip3 install specfabpy`
- Otherwise, for a local-only install, run `make specfabpy` in `/src` (requires LAPACK and BLAS)
- To install for general use (in other folders), run `make python` in `/src`. Note that if you do not have write permissions for your python installation, you can instead run `make specfabpy; python setup.py install --user`.

### Fortran
- The Fortran module is built by running `make specfab.o`

### Elmer/Ice Interface
- To interface with Elmer/Ice, you need a shared version of the libraries (built with the same Fortran compiler as Elmer). If needed, edit the compiler in the Makefile, then run `make libspecfab.so`.

## Q&A
- **Q** What *L* are possible?
  - **A** Any 4<=*L*<=20. If higher *L* are required:
    1. `cd src/include && python3 make_gaunt_coefs.py L` (replacing *L*)
    2. `make clean && make`

## Documentation
**Lattice rotation:** Rathmann et al. (2021), JOG, doi:10.1017/jog.2020.117 <br>
**Discontinuous dynamic recrystallization:** Rathmann and Lilien (2021), JOG, doi:10.1017/jog.2021.88 <br>
**Orthotropic bulk rheologies:** Rathmann and Lilien (in prep.)
