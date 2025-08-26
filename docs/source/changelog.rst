Changelog
=========

This page documents the changes made between releases of SALSA,
using the 1.0.0 release as the initial reference.

Version 2.0.0
-------------

* Overhauled the :class:`~salsa.AbsorberExtractor` API, including no longer requiring a ray filename to instantiate an object. See the :ref:`annotated-examples` for more.
* Switched :class:`~salsa.AbsorberExtractor` to a class inheritance framework so new methods will be easier to add as child classes (see :class:`~salsa.SPICEAbsorberExtractor`).
* Removed support for `spectacle`_ as it is no longer actively maintained and difficult to install.
* Overhauled the :class:`~salsa.AbsorberPlotter` API to match changes made to :class:`~salsa.AbsorberExtractor` and the removal of spectacle.
* Switched from using Pandas ``DataFrames`` to Astropy ``QTables`` in order to better track units.
* Removed ``cut_region`` functionality as it is no longer works with ``yt>=4.2``
* Stored the impact parameter of generated rays in the ``light_ray_solution`` list. See :ref:`using-lightray-solution` for more.
* Added minimal regression testing using `pytest`_.
* Updated documentation for consistent style.
* Removed unused functions.
* General maintenance and quality of life tweaks.

.. _spectacle: https://pypi.org/project/spectacle/
.. _pytest: https://pytest.org/