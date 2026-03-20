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

If you have the (minimal) dependencies already installed then you can use pip to install
SALSA: ::

  $ pip install astro-salsa

.. note::
  To install some of the dependencies, the GCC compiler needs to be installed.
  This is not a problem for most machines but can raise an error for some.

Development Version
-------------------

For the latest development version, follow these instructions if you have the
dependencies already installed. The easiest way to get all of the additional development
dependencies is to use the ``environment.yml`` file as detailed in the :ref:`conda-install` section.

With the dependencies installed, you next need to clone the `SALSA repository`_ 
then enter the main directory and use pip to
install the package: ::

  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
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
You may install the dependencies on your own, and you may have some of them installed
already (i.e., numpy, matplotlib) depending on your environment. To make installation
easier, we advise using conda and a conda environment to install.

.. _conda-install:

Conda Environment
-----------------

One of the easiest ways to make sure you have all the right dependencies is to
use a conda environment. In the repository there is an ``environment-minimal.yml`` file
that has all the packages necessary to run ``salsa``. The ``environment.yml`` file
additionally specifies the packages needed for development (e.g., pytest, sphinx). 

To create the environment
you first need to install `conda`_
and then run the following: ::

  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
  $ conda env create --file environment-minimal.yml
  $ conda activate salsa-env

.. _conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Now you should be able to painlessly install SALSA as described above!

.. note::
  This installs mpi4py using pip. This is useful if you already have an MPI library installed; see
  :ref:`install-mpi4py` for more details.

.. _manual_install:

Manually Installing Dependencies
---------------------------------

If you have a different, preferred method for installing dependencies you are more
than welcome to go with that route. See the ``environment-minimal.yml`` and/or ``setup.py``
in the `SALSA repository`_ for specifics about what
packages you need. Below is some advice/guides to installing a couple of the
trickier packages.

.. _install-yt:

Install yt
*************

yt can be installed in a few different ways but one of the easier ways is by
using conda. First
install `conda`_
then run: ::

  $ conda install -c http://use.yt/with_conda/ -c conda-forge yt

For full details about the different ways you can install yt, see
`yt's documentation`_.

.. _yt's documentation: https://yt-project.org/doc/

.. _install-trident:

Install Trident
****************

Trident can be installed via ``pip install trident``. The first time trident runs
though, it downloads an ionization table. It is recommended that you run trident
right after you pip install. This will also do some tests to make sure trident
is running properly. ::

  $ python
  >>> import trident; trident.verify()

For more details see `Trident's documentation`_

.. _Trident's documentation: https://trident.readthedocs.io/

.. _install-mpi4py:

Install mpi4py
**************

`mpi4py`_ is a package that
enables use of MPI parallelism with Python. SALSA uses this to split up light ray
creation and absorber extraction across multiple processors which becomes necessary
for large numbers of light rays.

.. _mpi4py: https://mpi4py.readthedocs.io/en/stable/index.html

mpi4py can be installed either using pip or conda. It's useful to install with pip
if you already have an MPI library installed, such as `OpenMPI`_ : ::

  $ pip install mpi4py

.. _OpenMPI: https://www.open-mpi.org/

If you want to use conda to install mpi4py, you need to be careful because there
may be problems if you have an MPI library already installed. Otherwise just: ::

  $ conda install mpi4py
