#! /usr/bin/env python
# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2024

import glob
from  distutils.core import setup

if len(glob.glob('specfabpy/*.so')) == 0:
    print('No compiled specfabpy found. Run `make specfabpy`')
else:
    setup(name='specfabpy',
          version='2024.8.17',
          author="Nicholas M. Rathmann and David A. Lilien",
          author_email="rathmann@nbi.ku.dk",
          description="specfab Python module",
          url="https://github.com/nicholasmr/specfab",
          packages=['.'],
          package_data={'': ['specfabpy/specfabpy.cpython*.so', 'specfabpy/*.py', 'specfabpy/fenics/*.py']},
    )
