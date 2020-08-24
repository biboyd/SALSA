from setuptools import setup, find_packages
import os

#base dependencies
#NOTE: gcc compiler is needed to install some packages
dependecies=['numpy', 'yt', 'trident', 'spectacle', 'matplotlib', 'pandas',
             'scipy', 'mpi4py']

setup(name="salsa",
      version="0.1.0",
      description = ("Synthetic absorber catalog generator from astrophysical simulations"),
      author="Brendan Boyd",
      author_email="boydbre1@msu.edu",
      license="BSD 3-Clause",
      keywords = ["simulation", "spectra", "astronomy", "astrophysics"],
      url="https://github.com/biboyd/SALSA",
      packages=find_packages(),
      install_requires=dependecies)
