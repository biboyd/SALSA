from setuptools import setup, find_packages
import os

#base dependencies
dependecies=['numpy', 'yt', 'trident', 'spectacle', 'matplotlib', 'pandas', 'scipy']

#if not building on rtd then add mpi4py
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if not on_rtd:
    dependecies += 'mpi4py'

setup(name="salsa",
      version="0.0.0",
      author="Brendan Boyd",
      packages=find_packages(),
      install_requires=dependecies)
