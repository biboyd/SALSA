.. _installation:

Installation
============

Installing SALSA
^^^^^^^^^^^^^^^^^

SALSA is built off of quite a few different packages and some of these are
best to install before installing SALSA. Go to :ref:`dependencies-install` for more details
on how best to do this.

Stable Version
--------------

SALSA has a few dependencies that are best installed 

The easiest way to install SALSA is with pip ::

  $ pip install astro-salsa

.. note::
  To install some of the dependencies, the GCC compiler needs to be installed.
  This is not a problem for most machines but can raise an error for some.

Development Version
-------------------

For the latest development version, you will first need to clone the `SALSA repository`_,
then enter the main directory. You can then install SALSA with pip: ::

  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
  $ pip install -e .

If you want to install the additional dependencies needed for locally testing SALSA
and building the documentation, instead run: ::

  $ pip install -e .[dev]

.. _SALSA repository: https://github.com/biboyd/SALSA

This will check that dependencies are installed and should install any that you
might be missing. The ``dev`` tag will check for
the additional packages needed to develop locally (e.g., pytest, sphinx).

Running the Test Suite
**********************

SALSA uses `pytest`_ for its test suite infrastructure.
To run the tests from either the top level ``SALSA`` directory or the ``tests`` subdirectory,
simply execute: ::

  $ pytest

.. _pytest: https://docs.pytest.org/

This test suite makes use of yt's `sample datasets`_.
By default this sample data will be downloaded to the directory from which the tests are run.
To keep things clean, you may wish to use the provided ``tests/data/`` directory
by setting yt's ``test_data_dir`` `configuration`_.
For example, you can locally configure ``test_data_dir`` within the SALSA directory using: ::

  $ yt config set --local yt test_data_dir ./tests/data/

Currently, the test suite is not written to utilize unit tests but rather it
regression tests SALSA's functionality wholistically.

.. _sample datasets: https://yt-project.org/doc/examining/loading_data.html#sample-data
.. _configuration: https://yt-project.org/doc/reference/configuration.html

.. _dependencies-install:

Installing Dependencies
^^^^^^^^^^^^^^^^^^^^^^^
For the full list of dependencies, see the ``requirements.txt`` and/or ``setup.py`` files
in the `SALSA repository`_. Below is some advice/guides to installing a couple of the
trickier packages.

A Note on Conda
---------------

`Conda <https://docs.conda.io/en/latest/>`_ is a popular package manager that
prefers to install wholly contained environments. This can cause issues with packages
like :ref:`mpi4py <install-mpi4py>` or h5py (often needed by :ref:`yt <install-yt>`)
that depend on non-Python libraries, as Conda
will install its own copy of these libraries. On most HPC systems, these libraries
are already provided and are usually optimized to the system itself.

It's generally discouraged to mix Conda and pip because it can cause issues with
dependency resolution, but it is required to install both SALSA and its dependency 
:ref:`Trident <install-trident>`. Additionally, 
of you wish to use a Conda environment but want to avoid duplicating non-Python libraries,
can often pip install such Python libraries using the existing external dependencies
in your system path.

Always make sure you are using the copy of pip associated with your Conda environment:

  $ conda activate my-env
  $ conda install pip
  $ which pip

.. _install-yt:

Install yt
----------

The yt package offers instructions on how to install it on 
`their website <https://yt-project.org/doc/installing.html#install-stable>`_. 

Note that if you want to work with a dataset from a particular simulation code,
there may be additional dependencies that need to be installed. For example, Enzo
datasets require you to install `h5py <https://docs.h5py.org/en/latest/build.html>`_.
This in turn depends on having `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`.

.. _install-trident:

Install Trident
---------------

Trident can be installed via pip: ::

  $ pip install trident
  
The first time trident runs
though, it downloads an ionization table. It is recommended that you run trident
right after you pip install. This will also do some tests to make sure trident
is running properly. ::

  $ python
  >>> import trident; trident.verify()

For more details see `Trident's documentation`_

.. _Trident's documentation: https://trident.readthedocs.io/

.. _install-mpi4py:

Install mpi4py
--------------

The `mpi4py`_ package enables use of MPI parallelism with Python. SALSA uses this to split up light ray
creation and absorber extraction across multiple processors which becomes necessary
for large numbers of light rays.

.. _mpi4py: https://mpi4py.readthedocs.io/en/stable/index.html

You can install mpi4py can be installed either using pip or Conda. It's useful to install with pip
if you already have an MPI library installed (such as `OpenMPI`_) as it will be 
built against your existing installation as long as the MPI compiler is in your system path: ::

  $ pip install mpi4py

.. _OpenMPI: https://www.open-mpi.org/

