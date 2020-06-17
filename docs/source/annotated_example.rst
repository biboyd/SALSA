.. _annotated-example:

Annotated Example
==================

To understand how salsa works and what it can do, here is an annotated example
to walk through a use case of salsa.

.. _catalog-generation:

Generate a Catalog
-------------------

One of the main goals of salsa is to make it easy to construct a catalog of
absorbers that can then be further analyzed. The process of constructing absorbers
takes two steps.

Step 1: Generate Light Rays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

trident is used to generate light rays that pass through your simulation. These
are 1-Dimensional objects that connect a series of cells and save field information
contained in those cells, ie density, temperature, line of sight velocity.
From these we can then extract absorbers along these light rays. (for further
information see `trident's documentation <https://trident.readthedocs.io/>`_)

To aid in generating these light rays, salsa contains the generate_lrays function
which can generate any number of lightrays which uniformly, randomly sample
impact parameters. This gives us a sample that is consistent with what the sample
of observational studies. This prevents any sampling bias when doing comparisons


Step 2: Extract Absorbers
^^^^^^^^^^^^^^^^^^^^^^^^^^

For more details on absorber extraction see :ref:`absorber-extraction`. But the
quick synopsis is there is the ice method and the spectacle method. Ice looks at
cell level data while spectacle fits lines to a synthetic spectra that is generated
by trident. The AbsorberExtractor class can use both methods.
