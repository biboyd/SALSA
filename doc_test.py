import numpy as np
import yt

import salsa

ds = yt.load_sample("HiresIsolatedGalaxy")

# define the center of the galaxy
center = [0.53, 0.53, 0.53]

# the directory where LightRays will be saved
ray_dir = 'my_rays'
n_rays = 4

# Choose what absorption lines to add to the dataset as well as additional
# field data to save
ion_list = ['H I', 'C IV']
other_fields = ['density', 'temperature', 'metallicity']
units_dict = {'density': 'g/cm**3',
              'temperature': 'K',
              'metallicity': 'Zsun'
              }

# the maximum distance a LightRay will be created (minimum default to 0)
max_impact = ds.quan(15, 'kpc')

np.random.seed(18)

salsa.generate_lrays(ds, center, n_rays, max_impact,
                    ion_list=ion_list, fields=other_fields, 
                    ray_directory=ray_dir)
#
# Plot code 
#
abs_ext = salsa.SPICEAbsorberExtractor(ds, ion_name='H I')
abs_ext.load_ray("my_rays/ray0.h5")
tab = abs_ext.get_current_absorbers(fields=other_fields, units_dict=units_dict)
print(tab)

ray_list = [f"{ray_dir}/ray0.h5",
            f"{ray_dir}/ray1.h5",
            f"{ray_dir}/ray2.h5",
            f"{ray_dir}/ray3.h5"]

# initialize a new AbsorberExtractor for looking at C IV
abs_ext_civ = salsa.SPICEAbsorberExtractor(ds, ion_name='C IV')
table = abs_ext_civ.get_all_absorbers(ray_list, fields=other_fields, units_dict=units_dict)
print(table)

# set the y limits for one of the plots
num_dense_min = 1e-11
num_dense_max = 1e-5

# Plot!
plotter = salsa.AbsorberPlotter(abs_ext)

fig, axes = plotter.plot_multiplot(outfname='example_multiplot.png',
                                   center = [0.53, 0.53, 0.53],
                                   num_dense_max=num_dense_max,
                                   num_dense_min=num_dense_min,
                                   make_spectra=False)