from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
#base dependencies
#NOTE: gcc compiler is needed to install some packages
dependecies=[
    'astropy',
    'h5py',
    'matplotlib',
    'mpi4py',
    'numpy',
    'scipy',
    'trident',
    'yt',
    ]

dev_deps=[
    'pooch',
    'pytest',
    'sphinx',
    'sphinx_automodapi',
    'sphinx_rtd_theme',
    ]

setup(name="astro-salsa",
      version="1.0.0",
      description = ("Synthetic absorber catalog generator from astrophysical simulations"),
      long_description=long_description,
      long_description_content_type='text/markdown',
      author="Brendan Boyd",
      author_email="boyd.brendan@stonybrook.edu",
      license="BSD 3-Clause",
      keywords = ["simulation", "spectra", "astronomy", "astrophysics"],
      url="https://github.com/biboyd/SALSA",
      packages=find_packages(),
      classifiers=["Programming Language :: Python :: 3.11",
                   "Programming Language :: Python :: 3.12",
                   "Programming Language :: Python :: 3.13"],
      python_requires=">=3.11, <3.14",
      install_requires=dependecies,
      extras_require={'dev':dev_deps})
