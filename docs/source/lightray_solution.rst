.. _using-lightray-solution:

Accessing Ray Spatial Information
=================================

It can be helpful to know where an extracted ray is located in relation to its parent dataset.
Thankfully, due to it's reliance on `Trident`_ for ray generation,
rays written to disk by SALSA contain information about their
generation. This information includes the start and end coordinates
of the ray as well as (unique to SALSA) the impact parameter used
to generate the rays. [#f1]_

.. _Trident: https://trident-project.org/

To access this information, we must first load a ray with `yt`_
Thanks to Trident and SALSA, this dataset will have a special attribute called ``light_ray_solution``.
This object is a list of dictionaries containing information about how the ray was generated. [#f2]_

.. _yt: https://yt-project.org/

.. code-block:: python

    import yt

    ds = yt.load("my_rays/ray0.h5")

    print(ds.light_ray_solution[0])

::

  {'end': unyt_array([0.52023164, 0.60173583, 0.55577773], 'unitary'),
   'filename': 'DD0044',
   'impact_parameter': unyt_quantity(0.00922918, 'code_length'),
   'redshift': 0.0,
   'start': unyt_array([0.53940663, 0.46652058, 0.48771738], 'unitary'),
   'traversal_box_fraction': unyt_quantity(0.15258786, 'unitary'),
   'unique_identifier': '1382182007'}

We can see that each key in the ``light_ray_solution`` dictionary has units attached.
These units (such as ``unitary`` and ``code_length``) may reference the coordinate system
of the simulation from which the light ray was extracted.

.. [#f1] Trident automatically writes start and end coordinates to its rays. 
    SALSA will add the impact parameter to any ray generated with :func:`~salsa.construct_rays`,
    which includes the two main ray generation mechanisms, :func:`~salsa.generate_lrays` and :func:`~salsa.generate_catalog`

.. [#f2] The ``light_ray_solution`` object is a list because Trident supports rays with multiple segments
    called `compound rays`_. Each segment can have it's own dictionary in ``light_ray_solution``.
    SALSA does not use these ray objects (instead only using "simple" rays) so
    ``light_ray_solution`` will only ever contain one dictionary.

.. _compound rays: https://trident.readthedocs.io/en/latest/annotated_example.html#compound-lightrays