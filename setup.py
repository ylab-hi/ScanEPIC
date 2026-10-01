#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build script for ScanEPIC's Cython extensions.

All package metadata (name, version, dependencies, entry points) lives in
pyproject.toml; this file only declares the compiled extensions.

@author: Josh Fry, Northwestern University, YangLab
"""

import numpy
import pysam
from Cython.Build import cythonize
from setuptools import Extension, setup

extensions = [
    Extension(
        name='src.extract.short._cython_fi',
        sources=['src/extract/short/_cython_fi.pyx'],
        extra_link_args=pysam.get_libraries(),
        include_dirs=pysam.get_include() + [numpy.get_include()],
        define_macros=pysam.get_defines()
    ),
    Extension(
        name='src.extract.single.cython_helpers',
        sources=['src/extract/single/cython_helpers.pyx'],
        extra_link_args=pysam.get_libraries(),
        include_dirs=pysam.get_include() + [numpy.get_include()],
        define_macros=pysam.get_defines()
    )
]

setup(
    ext_modules=cythonize(extensions, language_level=3),
)
