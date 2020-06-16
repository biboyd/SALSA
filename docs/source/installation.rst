.. _installation:

How To Install
==============


Step 1: Install YT
-------------------

The first step is to install YT. For full details, see
`YT documentation <https://yt-project.org/doc/>`_ .
One of the easiest way to install however is: ::

  $ conda install -c http://use.yt/with_conda/ -c conda-forge yt

Step 2: Install trident
------------------------

Trident can be installed via `pip install trident`. The first time trident runs
though, it downloads an ionization table. It is recommended that you run after
you pip install. ::

  $ python
  >> import trident; trident.verify()

For more details see `trident's documentation <https://trident.readthedocs.io/>`

Step 3: Install salsa
----------------------

Now you can install salsa. First `clone the repository <https://github.com/biboyd/SALSA>`_,
then enter the main directory and pip install: ::

  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
  $ pip install -e .

This should install any dependencies that you might be missing and allow you to
import the package.

Now you're all set to go! ::

  $ python
  >> import salsa
