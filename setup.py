#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name="spigen",
    version='0.1.0',
    author="Alexa Villaume",
    author_email="avillaum@ucsc.edu",
    packages=["spigen"],
    url="",
    license="LICENSE",
    description="Tools for generating stellar spectra from stellar parameters",
    #long_description=open("README.rst").read() + "\n\n"
                    #+ "Changelog\n"
                    #+ "---------\n\n"
                    #+ open("HISTORY.rst").read(),
    package_data={"spigen": ["Coefficients/*dat"]},
    include_package_data=True,
    #install_requires=["numpy", "scipy >= 0.9", "astropy", "matplotlib", "scikit-learn"],
)
