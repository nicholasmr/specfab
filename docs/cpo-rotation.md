# Rotation

Rotating an expansion series by `theta` about the $y$-axis (in the $x$&mdash;$z$ plane) followed by `phi` about the $z$-axis (in the $x$&mdash;$y$ plane) is done by:

```python
import numpy as np
from specfabpy import specfabpy as sf
lm, nlm_len = sf.init(8) 

### Construct an arbitrary fabric to rotate
a2 = np.diag([0.0,0.0,1.0]) # any second-order structure tensor (not necessarily diagonal)
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients of corresponding normalized ODF

### Rotate ODF
# Assumes L=<12 (rotation for larger L is not implemented)
theta = np.deg2rad(-45) 
phi   = np.deg2rad(45)
nlm_rot1 = sf.rotate_nlm(nlm, theta, 0)    # first rotate around y axis in x-z plane
nlm_rot2 = sf.rotate_nlm(nlm_rot1, 0, phi) # next  rotate around z axis in x-y plane 
```

Plotting the expansion series (`nlm` and `nlm_rot*`) gives:

![](https://github.com/nicholasmr/specfab/raw/main/tests/wigner-d-rotation-test/wigner-d-rotation-test.png){: style="width:600px"}
