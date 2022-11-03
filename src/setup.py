#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 dlilien <dlilien@hozideh>
#
# Distributed under terms of the MIT license.

"""
Install the pre-compiled so for python
"""
import glob
from distutils.core import setup

if len(glob.glob('specfabpy.cpython*.so')) == 0:
    print('No compiled specfabpy found. Run `make python`')
else:
    setup(name='specfabpy',
          author="Nicholas M. Rathmann and David A. Lilien",
          author_email="rathmann@nbi.ku.dk",
          description="specfab Python module",
          url="https://github.com/nicholasmr/specfab",
          version='2022.11.3',
          packages=['.'],
          package_data={'': ['specfabpy.cpython*.so']},
          )
