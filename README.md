## specfab

<img src="https://raw.githubusercontent.com/nicholasmr/specfab/main/images/logo-square.jpg" alt="specfab logo" width="330px"> 

Spectral CPO model of polycrystalline materials that:

- Can model lattice rotation, discontinuous DRX, and rotation/continuous DRX.
- Can calculate CPO-induced viscous anisotropies using Sachs/Taylor homogenizations.
- Can calculate elastic P- and S-wave velocities using Voigt/Reuss homogenizations.
- Contains expressions for forward+inverse orthotropic and transversely isotropic rheologies.
- Can convert between structure tensors and spectral expansions coefficients.
- Can be integrated with finite-element models such as Elmer and FEniCS.

## Documentation

See [the specfab docs](https://nicholasmr.github.io/specfab) for installation, tutorials, and more.

## Q&A

| **Q** | **A** |
| :--- | :--- |
| What $L$ are possible? | Any $4\leq L\leq 20$. If higher $L$ are required: <br>1. `cd src/include && python3 make_gaunt_coefs.py L` (replacing `L`) <br>2. `make clean && make` |
