# SALSA
[![Build Status](https://travis-ci.com/biboyd/SALSA.svg?branch=master)](https://travis-ci.com/biboyd/SALSA)
[![Documentation Status](https://readthedocs.org/projects/salsa/badge/?version=latest)](https://salsa.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/biboyd/SALSA/master?filepath=notebooks%2FExample_notebook.ipynb)

SALSA: Synthetic Absorption Line Surveyor Application. Python tool that combines YT, Trident, Spectacle and other python packages to create a steamlined pipeline to construct absorber catalogs from galactic simulations. Currently only ENZO datasets have been used/tested.

Read the Docs [here](https://salsa.readthedocs.io) 

## Install
First clone the repo, enter the directory and run 
```
  $ git clone https://github.com/biboyd/SALSA.git
  $ cd SALSA
  $ pip install -e .
  $ python
  >>> import salsa
``` 
Now you should be all set to code!  
For more detailed description see the [installation guide](https://salsa.readthedocs.io/en/latest/installation.html). 

## Getting Started
For an annotated example [go here](https://salsa.readthedocs.io/en/latest/annotated_example.html). Or lanuch an interactive jupyter hosted on Binder, [here](https://mybinder.org/v2/gh/biboyd/SALSA/master?filepath=notebooks%2FExample_notebook.ipynb).

The easiest way to get started is use `salsa.generate_catalog()`. This simply takes the dataset, number of light rays to make, directory to save those light rays, a list of ions and some other optional parameters. This creates a number lightrays and then extracts absorbers for each ion and returns a `pandas.DataFrame` that can then be further analyzed.

