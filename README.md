# Spectral orientation fabric 

Model by Rathmann et al. (2020)

**Contains:**
- Modules for: Fortran (`make specfab.o`) and Python (`make specfabpy`).
- Demos for: Fortran (`make demo`) and Python (`make demopy`).
- `plot.py` script for plotting solutions of the demos (netCDF dumps).

## Q&A
- **Q** What *L* are possible?
  - **A** 8<=*L*<=40 are allowed. If higher *L* are required:
    1. `cd include && python3 make_gaunt_coefs.py L` (replacing *L*)
    2. Edit `gaunt.f90` accordingly
    3. `make clean && make specfab.o`
