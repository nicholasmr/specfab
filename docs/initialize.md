# Initialization

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

