.. _parallelism:

Running SALSA in Parallel
=========================

To get a statistically meaningful sample, SALSA should be used to generate
a large number of rays. Extracting rays from a dataset and then finding
absorbers within these rays can be a computationally expensive process.
Thankfully, since each ray is independent of each other ray, this process
can be easily parallelized.

SALSA is already equipped to generate rays in parallel using `MPI`_.
Additionally, :func:`~salsa.generate_catalog` can be used
to extract absorbers in parallel. **This means that running SALSA in parallel
is as simple as having multiple CPUs available to run on!**

.. _MPI: https://en.wikipedia.org/wiki/Message_Passing_Interface

.. warning::
    It's not advisable to use SALSA within a Python script that itself uses
    MPI.

The following functions are written to be able to run in parallel:

* :func:`salsa.generate_lrays`
* :func:`salsa.construct_rays`
* :func:`salsa.generate_catalog` (both ray generation and absorber extraction)
* :func:`salsa.utils.check_rays`