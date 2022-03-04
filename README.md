# Spectral orientation fabric model for polycrystals
A spectral fabric model for polycrystalline materials (e.g. ice, olivine, etc.) that 
- Can model fabric processes such as kinematic lattice (c-axis) rotation (Svendsen and
Hutter, 1996), discontinuous dynamic recrystallization (Placidi and others, 2010), and rotation recrystallization (Gödert, 2003).
- Can calculate fabric-induced directional viscosities (directional strain-rate enhancement factors) using Sachs/Taylor homogenization schemes.
- Contains analytical expressions for both forward and inverse anisotropic rheologies such as the transversely isotropic rheology and the orthotropic rheology.

## Demo of lattice rotation and resulting strain-rate enhancement factors
![image](demo/cube-crush-animation/cube-crush.gif)

## Installation
### Python module `specfabpy`
- If you are running Linux, a pre-compiled version is available: `pip3 install numpy --upgrade && pip3 install specfabpy`
- Else run `make specfabpy` in `/src` (requires LAPACK and BLAS)
### Fortran
- Fortran module is built by running `make specfab.o`

## Demos
- Basic demos using `specfabpy` can be found on the [Wiki](https://github.com/nicholasmr/specfab/wiki)
- More demos are located in `/demo`. Note that Fortran demos must be built seperately (see `Makefile`), e.g. `make fabric-evolution-demo` .

## Q&A
- **Q** What *L* are possible?
  - **A** Any 4<=*L*<=20. If higher *L* are required:
    1. `cd src/include && python3 make_gaunt_coefs.py L` (replacing *L*)
    2. `make clean && make`

## Documentation
**Lattice rotation:** Rathmann et al. (2021), JOG, doi:10.1017/jog.2020.117 <br>
**Discontinuous dynamic recrystallization:** Rathmann and Lilien (2021), JOG, doi:10.1017/jog.2021.88 <br>
**Orthotropic bulk rheologies:** Rathmann and Lilien (in prep.)
