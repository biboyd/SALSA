.. SALSA documentation master file, created by
   sphinx-quickstart on Tue Jun 16 15:17:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SALSA's documentation!
=================================

The Synthetic Absorption Line Surveyor Application (SALSA)
is an open-source Python module that creates a streamlined process to
generate synthetic absorber catalogs from galactic simulations. Multiple open-source software
projects utilized to achieve this: 
accessing simulation data is done using `yt`_,
while `Trident`_ is used to generate synthetic 
sight lines/light rays and synthetic spectra. 
`Astropy`_ ``QTables`` are used to store data about identified absorbers.

.. _yt: https://yt-project.org/
.. _Trident: https://trident-project.org/
.. _Astropy: https://www.astropy.org/

Observational studies generate large absorber catalogs by studying the absorption
line spectra of distant quasars.
As quasar light passes through intervening galaxies,
each galaxy leaves its "imprint" on the quasar spectrum by absorbing some wavelength of light.
SALSA can generate similar catalogs from cosmological and galactic simulations,
allowing research to study these simulations from an observer's perspective. This
can give new insights into the data as well as help facilitate comparisons and
collaboration between simulations and observations

SALSA includes a novel method for extracting absorbers, the SPICE method. This uses
cell-level data to extract absorbers from a Trident ``LightRay`` and returns a great
deal of information that can be further analyzed.

.. toctree::
   :maxdepth: 2

   installation
   annotated_example
   absorber_extraction
   parallelism
   lightray_solution
   reference
   changelog

Index
==================
* :ref:`genindex`
