Changelog
=========

This page documents the changes made between releases of SALSA,
using the 1.0.0 release as the initial reference.

Version ???
-----------

* Removed support for `spectacle <https://pypi.org/project/spectacle/>`_ as it is no longer actively maintained and difficult to install.
* Switched from using Pandas ``DataFrames`` to Astropy ``QTables`` in order to better track units.
* Added minimal regression testing using `pytest <pytest.org>`_
* Updated documentation for consistent style
* General maintenance to keep up with package dependencies, etc.