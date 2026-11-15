"""Compatibility entry point for setuptools.

Package metadata, dependencies, and the cpu/cuda/rocm/docs installation extras
are declared in pyproject.toml so all installation methods use the same values.
"""

from setuptools import setup


setup()
