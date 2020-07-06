.. SALSA documentation master file, created by
   sphinx-quickstart on Tue Jun 16 15:17:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SALSA's documentation!
=================================

Salsa is an open-source python module that creates a streamlined process to
generate synthetic absorber catalogs from galactic simulations. Multiple open-source software
projects utilized to achieve this. Accessing simulation data is done using yt.
Trident is used to generate synthetic sightlines/lightrays and generate synthetic
spectra. Spectacle is used to fit voigt profiles to spectra and extract absorbers.

In addition a novel method for extracting absorbers, the Ice method. This uses
cell level data to extract absorbers from a Trident lightray and returns a great
deal of information that can be further analyzed.

.. toctree::
   :maxdepth: 2

   installation.rst
   annotated_example.rst
   absorber_extraction.rst
   reference.rst

Index
==================
* :ref:`genindex`
