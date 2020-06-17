.. _annotated-example:

Annotated Example
==================

To understand how salsa works and what it can do, here is an annotated example
to walk through a use case of salsa.

.. _extract-absorbers-example:

Extracting Absorbers
---------------------

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

To aid in generating these light rays, salsa contains the
:class:`~salsa.generate_lrays` function which can generate any number of lightrays
which uniformly, randomly sample impact parameters. This gives us a sample that
is consistent with what the sample of observational studies. This prevents any
sampling bias when doing comparisons.

To use this function:::

  import yt
  import salsa
  import numpy as np
  import pandas as pd

  # load in the simulation dataset
  ds_file = "IsolatedGalaxy/galaxy0030/galaxy0030"
  ds = yt.load(ds_file)

  # define the center of the galaxy
  center= ds.arr([0.5, 0.5, 0.5], 'code_length')

  # the directory where lightrays will be saved
  ray_dir = 'my_rays'
  n_rays=4

  # Choose what absorption lines to add to the dataset as well as additional
  # field data to save
  ion_list=['H I', 'C IV', 'O VI']
  other_fields = ['density', 'temperature', 'radius']

  # the maximum distance a lightray will be created (minimum default to 0)
  max_impact = 200 #kpc

  # set a seed so the function produces the same random rays
  np.random.seed(42)

  # Run the function and rays will be saved to my_rays directory
  salsa.generate_lrays(ds, center, n_rays, max_impact, center=center,
                       ion_list=ion_list, fields=other_fields, out_dir=ray_dir)


Step 2: Extract Absorbers
^^^^^^^^^^^^^^^^^^^^^^^^^^

For more details on absorber extraction see :ref:`absorber-extraction`. But the
quick synopsis is there is the ice method and the spectacle method. Ice looks at
cell level data while spectacle fits lines to a synthetic spectra that is generated
by trident. The :class:`~salsa.AbsorberExtractor` class can use both methods.

Now let's extract some absorbers from the Light rays we made
::

  ray_file = f"{ray_dir}/ray0.h5"

  # construct absorber extractor
  abs_ext = salsa.AbsorberExtractor(ds, ray_file, ion_name='H I')

  # use ice method to extract absorbers into a pandas DataFrame
  units_dict=dict(density='g/cm**3', radius='kpc')
  df_ice = abs_ext.get_ice_absorbers(other_fields, user_unit_dict=units_dict)
  df_ice.head()

  # use spectacle now
  df_spect = abs_ext.get_spectacle_absorbers()
  df_spect.head()

Notice that both of these methods contain different information. Ice includes
more details of the simulation data like the density and temperature of the
absorber, something that is not easily detected from the spectra. Spectacle
contains more information of the line like the equivalent width and the doppler
b parameter.

To extract absorbers from multiple ``LightRays`` you can use the
:class:`~salsa.get_absorbers` function. This will loop through a list of rays and
extract absorbers from each one. see:::

  ray_list = [f"{ray_dir}/ray0.h5",
              f"{ray_dir}/ray1.h5",
              f"{ray_dir}/ray2.h5",
              f"{ray_dir}/ray3.h5"]

  # initialize a new AbsorberExtractor for looking at O VI
  abs_ext_ovi = salsa.AbsorberExtractor(ds, ray_file, ion_name='O VI')

  df_ovi = get_absorbers(abs_ext_ovi, ray_list, method='ice',
                         fields=other_fields, user_unit_dict=units_dict)

  df_ovi.head()

Notice that the spectacle method could also be used. Also, although the
AbsorberExtractor takes a ray file at construction, new rays can be loaded into
it.

To retain information on where each absorber came from, an ``absorber_index`` is
given. The number represents the ray it was extracted from and the letter
signifies the order in which the absorber was extracted. So the first absorber
to be extracted from ray2.h5 would have an index of ``2A`` and the next would be
``2B``. This can be useful for comparing/analyzing absorbers on the same sightline.

.. _catalog-generation-example:

Catalog Generation
-------------------
To generate a full catalog of absorbers we can use the
`:class:~salsa.generate_catalog` function to both generate a sample of
``trident.LightRay`` objects and then :class:`~salsa.AbsorberExtractor` to extract
absorbers of a list of ions.

Here is what you need to setup and run:::

  df_catalog = salsa.generate_catalog(ds, n_rays, ray_dir, ion_list,
                                      fields=other_fields, center=center,
                                      impact_param_lims=(0, max_impact),
                                      method='ice', units_dict=units_dict)

  df_catalog.head()

This function looks first to see if rays have been created in the given directory.
If there are the right number of rays and they all contain the right ions and
other fields that were specified (in this case that would be 'density',
'temperature', 'radius'), then those rays will be used. Otherwise, new rays are
created using :class:`~salsa.generate_lrays`.

Next, :class:`~salsa.get_absorbers` is used to find the absorbers from each ion
in ``ion_list`` and finally a catalog is returned as a ``pandas.DataFrame``. Note
that the absorber index is unique only up to the ion/wavelength


.. _visualizing-absorbers:

Visualizing Absorbers
---------------------
To visualize what is actually be extracted from the ``LightRay`` objects and
synthetic spectra, you can use the :class:`~salsa.AbsorberPlotter` class. This
is built off of the :class:`~salsa.AbsorberExtractor` with added functionality
to make plots.

To get a full picture of what is happening at each level we can create a
multi panel plot containing:

    1. a slice of the simulation with the ray annotated
    2. The number density profile along the ray's path
    3. The line of sight velocity profile along the ray's path
    4. The synthetic spectra created from the ray

This figure gives you a good overview of what is happening and can give valuable
context to the absorption extraction methods. Additionally, each plot can be made
individually if you care less about the spectra, or don't want to plot a slice
(which can be time consuming, depending on the detail in the simulation).

To create the multi-panel plot:::

  import salsa
  import yt
  import matplotlib.pyplot as plt

  ds_file="IsolatedGalaxy/galaxy0030/galaxy0030"

  # using one of the rays generated in previous example.
  # Any light ray can be used.
  ray_file = "my_rays/ray0.h5"

  plotter = salsa.AbsorberPlotter(ds_file, ray_file, "H I",
                                  center_gal=[0.5, 0.5, 0.5],
                                  use_spectacle=True,
                                  plot_spectacle=True,
                                  plot_ice=True)

  fig, axes = plotter.create_multi_plot()

  fig.show()

The grey regions on the middle two plots indicate the absorbers that the ice
method finds. The three highest column densities are marked and displayed in a
legend. In the last plot, the solid lines indicate the "raw" spectra while the
dotted lines show the absorption lines that spectacle fit (only the three largest
lines are plotted with their column densities recorded in a legend).

The total column density along the lightray, the total found via the ice method
and the total found by spectacle is recorded in a legend in the spectra plot.
