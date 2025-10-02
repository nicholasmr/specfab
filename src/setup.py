#! /usr/bin/env python
# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-

import glob
from setuptools import setup, find_packages

if len(glob.glob('specfabpy/*.so')) == 0:
    print('No compiled specfabpy found. Run `make specfabpy`')
else:
    setup(
        name="specfabpy",
        version="2025.10.2",
        author="Nicholas M. Rathmann and David A. Lilien",
        author_email="rathmann@nbi.ku.dk",
        description="specfab Python module",
        url="https://github.com/nicholasmr/specfab",
        install_requires=[
            "numpy>=1.24", #"scipy", "matplotlib", "cmasher", "cartopy",
        ],
        packages=find_packages(include=["specfabpy", "specfabpy.*"]),
        include_package_data=True,
        package_data={
            "specfabpy": [
                "*.so", "*.py",
                "fenics/*.py",
                "firedrake/*.py",
                "tamm/*.py",
            ],
        },
        python_requires=">=3.8",
    )

