#!/usr/bin/env python
# All metadata lives in setup.cfg; the version comes from the git tags (setuptools_scm).
from setuptools import setup

setup(use_scm_version={"fallback_version": "0.0.0"})
