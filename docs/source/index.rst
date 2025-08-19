.. SALSA documentation master file, created by
   sphinx-quickstart on Tue Jun 16 15:17:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SALSA's documentation!
=================================

SALSA is an open-source Python module that creates a streamlined process to
generate synthetic absorber catalogs from galactic simulations. Multiple open-source software
projects utilized to achieve this: 
accessing simulation data is done using `yt <https://yt-project.org/>`_,
while `Trident <https://trident-project.org/>`_ is used to generate synthetic 
sight lines/light rays and synthetic spectra. 
Astropy <https://www.astropy.org/>`_ ``Tables`` are used to store data about identified absorbers.

Observational studies generate large absorber catalogs by studying the absorption
line spectra of distant quasars.
As quasar light passes through intervening galaxies,
each galaxy leaves its "imprint" on the quasar spectrum by absorbing some wavelength of light.
SALSA can generate similar catalogs from cosmological and galactic simulations,
allowing research to study these simulations from an observer's perspective. This
can give new insights into the data as well as help facilitate comparisons and
collaboration between simulations and observations

SALSA includes a novel method for extracting absorbers, the SPICE method. This uses
cell-level data to extract absorbers from a Trident ``lightray`` and returns a great
deal of information that can be further analyzed.

.. toctree::
   :maxdepth: 2

   installation.rst
   annotated_example.rst
   absorber_extraction.rst
   reference.rst
   changelog.rst

Index
==================
* :ref:`genindex`
