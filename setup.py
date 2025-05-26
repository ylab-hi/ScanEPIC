#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Setup script for ScanEPIC - A tool for identifying exitron splicing events from RNA-seq data

@author: Josh Fry, Northwestern University, YangLab
"""

from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize
import pysam
import numpy

name = "scanepic"
__version__ = "2.0.0"

# Define Cython extensions
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

# Read long description from README if it exists
try:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "ScanEPIC - A tool for identifying exitron splicing events from RNA-seq data"

setup(
    name=name,
    version=__version__,
    author="Josh Fry",
    author_email="your.email@northwestern.edu",  # Update with actual email
    description="A tool for identifying exitron splicing events from RNA-seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/scanepic",  # Update with actual repo
    packages=find_packages(),
    ext_modules=cythonize(extensions),
    entry_points={
        'console_scripts': [
            'scanepic=src.cli:cli',
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",  # Update as needed
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pysam>=0.19.0",
        "click>=8.0.0",
        "gffutils>=0.10.0",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "biopython>=1.79",
        "statsmodels>=0.12.0",
    ],
    extras_require={
        "dev": [
            "pytest",
            "black",
            "flake8",
        ],
        "docs": [
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
    include_package_data=True,
    package_data={
        'src.extract.long': ['blacklist.tsv'],
    },
    zip_safe=False,  # Important for Cython extensions
)
