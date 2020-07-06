.. _installation:

How To Install
==============

Installing Salsa
^^^^^^^^^^^^^^^^^

Salsa is built off of quite a few different packages so it is best to install
these before installing salsa. Go to :ref:`dependencies-install` for more details
on how best to do this.

If you have the dependencies already installed, you need to `clone the
repository <https://github.com/biboyd/SALSA>`_ if you haven't already then enter
the main directory and use pip to install the package: ::

  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
  $ pip install -e .

This will check that dependencies are installed and should install any that you
might be missing. Again, since this is a somewhat complicated environment, it is
best to install the dependencies beforehand.

Now you're all set to go! ::

  $ python
  >>> import salsa

.. _dependencies-install:

Installing Dependencies
^^^^^^^^^^^^^^^^^^^^^^^
You may install the dependencies on your own, and you may have some of them installed
already (ie numpy, matplotlib) depending on your environment. To make installation
easier, we advise using conda and a conda environment to install.

.. _conda-install:

Conda Environment
-----------------

One of the easiest ways to make sure you have all the right dependencies is to
use a conda environment. In the repository there is an ``enivronment.yml`` file
that has all the packages necessary to run ``salsa``. To create the enivronment
you first need to
`install conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
and then run the following: ::

  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
  $ conda env create --file environment.yml
  $ conda activate salsa-env

Now you should be able to painlessly install salsa as described above!

.. note::
  This installs mpi4py using conda. This may cause problems if you already have
  an MPI Library already installed because conda will try to install one itself.
  If you already have an MPI library it is best to use pip see
  :ref:`install-mpi4py` for more details.

.. _manual_install:

Manually Installing Dependencies
---------------------------------

If you have a different, preferred method for installing dependencies you are more
than welcome to go with that route. See the ``environment.yml`` and/or ``setup.py``
in the `Github repo <https://github.com/biboyd/SALSA>`_ for specifics about what
packages you need. Below is some advice/guides to installing a couple of the
trickier packages.

.. _install-yt:

Install yt
*************

yt can be installed in a few different ways but one of the easier ways is by
using conda. First
`install conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
then run: ::

  $ conda install -c http://use.yt/with_conda/ -c conda-forge yt

For full details about the different ways you can install yt, see
`yt's documentation <https://yt-project.org/doc/>`_ .

.. _install-trident:

Install Trident
****************

Trident can be installed via ``pip install trident``. The first time trident runs
though, it downloads an ionization table. It is recommended that you run trident
right after you pip install. This will also do some tests to make sure trident
is running properly. ::

  $ python
  >>> import trident; trident.verify()

For more details see `trident's documentation <https://trident.readthedocs.io/>`_

.. _install-mpi4py:

Install mpi4py
**************

`mpi4py <https://mpi4py.readthedocs.io/en/stable/index.html>`_ is a package that
enables use of MPI parallelism with python. Salsa uses this to split up lightray
creation and absorber extraction across multiple processors which becomes necessary
for large numbers of lightrays.

mpi4py can be installed either using pip or conda. If you want to pip install
mpi4py, you need to have an MPI library already installed, like
`OpenMPI <https://www.open-mpi.org/>`_ . Otherwise just: ::

  $ pip install mpi4py

If you want to use conda to install mpi4py, you need to be careful because there
may be problems if you have an MPI library already installed. Otherwise just: ::

  $ conda install mpi4py
