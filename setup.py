from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
#base dependencies
#NOTE: gcc compiler is needed to install some packages
dependecies=['trident~=1.3',
             'yt~=4.2',
             'matplotlib',
             'numpy',
             'scipy',
             'pandas',
             'mpi4py']

opt_spectacle=['astropy<=4',
               'specutils<=0.6.1', 
               'gwcs<=0.11.0', 
               'spectacle']

opt_testing=["pooch",
             "pytest"]

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
      classifiers=["Programming Language :: Python :: 3",
                   "Programming Language :: Python :: 3.8",
                   "Programming Language :: Python :: 3.9",
                   "Programming Language :: Python :: 3.10",
                   "Programming Language :: Python :: 3.11"],
      python_requires=">=3.8, <3.12",
      install_requires=dependecies,
      extras_require={'spectacle':opt_spectacle,
                      'testing':opt_testing})
