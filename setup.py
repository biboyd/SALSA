from setuptools import setup, find_packages
import os

#base dependencies
dependecies=['numpy', 'yt', 'trident', 'spectacle', 'matplotlib', 'pandas', 'scipy', 'mpi4py']

setup(name="salsa",
      version="0.0.0",
      author="Brendan Boyd",
      packages=find_packages(),
      install_requires=dependecies)
