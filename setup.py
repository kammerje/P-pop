#!/usr/bin/env python
# Minimal setup.py; most config information is now in pyproject.toml

from setuptools import setup, find_packages

setup(packages=find_packages(),
      use_scm_version=True,
      setup_requires=['setuptools_scm'],
      )
