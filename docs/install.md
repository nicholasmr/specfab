# Install

* Requires a version of Python between **3.6 and 3.11** &mdash; newer versions are not yet supported because specfab relies on `setuptools<=61`.
* Source code is [available here](https://github.com/nicholasmr/specfab).

## 🐧 Linux

The easiest way to install is using the *pipy* package:
```bash
pip3 install numpy --upgrade && pip3 install specfabpy
```
Alternatively, you can compile specfab yourself:
```bash
cd src
pip3 install setuptools==61
make install
```

## 🍎 Mac

You must compile specfab yourself:
```bash
cd src
pip3 install setuptools==61
make all
pip3 install . # or run: make install
```

## 🪟 Windows

No reported experience installing specfab on Windows. You will have to try this yourself, but please report back if successful, and how 🚀.

## 🐍 Conda

Conda users have reported this to work:
```bash
brew install gcc # if gcc is needed on Mac, else skip

conda create --name sfenv
conda install -c anaconda ipykernel # for jupyter support
python -m ipykernel install --user --name=sfenv
conda activate sfenv

conda install python==3.11
pip install numpy==1.23
pip install setuptools==61
pip install cartopy, scipy, cmasher

cd specfabpy/src
make install
```

## ⚙️ Fortran interfaces

| Interface | Make |
| :--- | :--- |
| Plain Fortran interface | Run `cd src && make specfab.o` |
| Elmer interface | Compile shared library by running `cd src && make libspecfab.so` |

*Libraries required:* BLAS, LAPACK<br>

# Initialize 

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

