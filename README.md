# SALSA
[![Build Status](https://github.com/biboyd/SALSA/actions/workflows/python-app.yml/badge.svg?branch=master&event=push)](https://github.com/biboyd/SALSA/actions/workflows/python-app.yml)
[![Documentation Status](https://readthedocs.org/projects/salsa/badge/?version=latest)](https://salsa.readthedocs.io/en/latest/?badge=latest)
<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/biboyd/SALSA/master?filepath=notebooks%2FExample_notebook.ipynb) -->
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02581/status.svg)](https://doi.org/10.21105/joss.02581)
[![DOI](https://zenodo.org/badge/271633933.svg)](https://zenodo.org/badge/latestdoi/271633933)

SALSA (Synthetic Absorption Line Surveyor Application) is a Python tool that
constructs synthetic absorber catalogs from hydrodynamic galaxy simulations.
SALSA heavily utilizes [yt](https://yt-project.org/) to access simulation data
and [Trident](http://trident-project.org/) to create light rays/sight lines and
generate synthetic spectra.

Observational studies generate large absorber catalogs by studying the absorption
line spectra of distant quasars, as their light passes through intervening galaxies.
SALSA can generate similar catalogs from cosmological and galactic simulations,
allowing research to study these simulations from an observers perspective. This
can give new insights into the data as well as help facilitate comparisons and
collaboration between simulations and observations.

SALSA allows us to dip into galactic simulations and start to chip away at the
many unknowns of the universe

A [JOSS paper](https://joss.theoj.org/papers/10.21105/joss.02581) was published for 
SALSA and we recommend reading it for an overview of the package and its possible uses. 
If you do use SALSA in a project we ask that you cite this paper.

For detailed information on how to install and run SALSA, Read the Docs
[here](https://salsa.readthedocs.io)

## Install

We recommending checking the [documentation](https://salsa.readthedocs.io/en/latest/installation.html#installing-dependencies>) for advice on installing the
core dependencies.

Once you have the trickier dependencies installed, you can use pip and
run these commands to install the most stable version:
```
  $ pip install astro-salsa
  $ python
  >>> import salsa
```

If you want to install the latest development version, you can clone the repository and run these commands to install
the additional dependencies necessary for testing and building
the documentation:
```
  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
  $ pip install -e .[dev]
  $ python
  >>> import salsa
```
Now you should be all set to code!

### Installing dependencies
The `requirements.txt` file contains the minimum dependencies needed to run SALSA.

Note that you need gcc compiler installed (which it often already is on most machines).
For a more detailed description see the
[installation guide](https://salsa.readthedocs.io/en/latest/installation.html)
which also includes tips if you want to install dependencies on your own.

To install the basic dependencies with pip, run

```
pip install -r requirements.txt
```

If you prefer to use conda, you cannot use the `requirements.txt` file
directly as [Trident](http://trident-project.org/) is not available via conda.
You'll also still need to install SALSA with pip as above. **Make sure `pip` is
installed in your conda environment!**

## Getting Started
For an annotated example [go here](https://salsa.readthedocs.io/en/latest/annotated_example.html). Or launch an interactive jupyter hosted on Binder
[here](https://mybinder.org/v2/gh/biboyd/SALSA/master?filepath=notebooks%2FExample_notebook.ipynb) (note that the notebook
may take some time to load as it generally has to build the repository).

If you want to explore on your own, the easiest way to get started is use
`salsa.generate_catalog()`. This takes:
  * The simulation dataset (either as a path or an object loaded by yt)
  * Number of light rays/sightlines to make
  * Directory to save those light rays
  * A list of ions
  * Some other optional parameters.  
This creates a number light rays and then extracts absorbers for each ion. An
Astropy [`QTable`](https://docs.astropy.org/en/stable/table/) is returned with information about all the absorbers which
can then be further analyzed.

### Building the Documentation
If you installed SALSA with
`pip install astro-salsa[dev]` or `pip install -e .[dev]`, you'll be able to build the documentation locally:
```
  $ cd docs
  $ make html  # or latexpdf if you prefer
```

This will create a `docs/build/` folder that contains the documentation.

## Contributing Guidelines
All contributions are welcome! This is an open-source project, built on many
other open-source projects. Contributing can take many forms including:
contributing code, testing and experimenting, or offering ideas for different
features.

If you are interested in contributing you can contact us directly at
kopenhaf@msu.edu or add an issue on this Github page.
