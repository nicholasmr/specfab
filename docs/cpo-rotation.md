# Rotation

Rotating an expansion series by $\theta$ about the $y$-axis (in the $x$&mdash;$z$ plane) followed by $\phi$ about the $z$-axis (in the $x$&mdash;$y$ plane) can be done as follow:

```python
import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(8) 

### Construct an arbitrary fabric to rotate
a2 = np.diag([0, 0, 1]) # arbitrary second-order structure tensor, a^(2)
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients of corresponding ODF

### Rotate ODF
# Note: assumes L=<12 (rotation for larger L is not implemented)
theta = np.deg2rad(-45) 
phi   = np.deg2rad(45)
nlm_rot1 = sf.rotate_nlm(nlm, theta, 0)    # first rotate around y axis in x-z plane
nlm_rot2 = sf.rotate_nlm(nlm_rot1, 0, phi) # next  rotate around z axis in x-y plane 
nlm_rot3 = sf.rotate_nlm(nlm_rot2, -theta, -phi) # rotate back

# See "plotting" pages on how to plot the resulting ODFs
```

![](https://github.com/nicholasmr/specfab/raw/main/tests/rotate-Wigner-D/wigner-d-rotation-test.png){: style="width:650px"}
