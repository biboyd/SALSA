import yt
import numpy as np
from sys import argv

#
#Find the center and orientation of galaxy
#return center coordinates and normal vector
#

def find_center(ds_fname, center=None):
    ds = yt.load(ds_fname)

    if center is None:
        #let max density be center of gal
        v, center = ds.find_max('density')
    else:
        center = ds.arr(center, 'code_length')
    #create sphere around disk of galaxy
    sph_gal = ds.sphere(center, (20, 'kpc'))

    #compute the bulk velocity
    bulk_velocity = sph_gal.quantities.bulk_velocity(use_particles=False)
    sph_gal.set_field_parameter(("bulk_velocity"), bulk_velocity)
    #compute angular momentum vector normalize it
    axis_rot = sph_gal.quantities.angular_momentum_vector()
    norm_vector = axis_rot/(np.linalg.norm(axis_rot) * axis_rot.unit_array)

    return center.in_units('code_length'), norm_vector


if __name__ == '__main__':
    ds_filename = argv[1]
    c, n = find_center(ds_filename)

    print("x, y, z")
    print(f"{c[0]}, {c[1]}, {c[2]}")
    print(f"{n[0]}, {n[1]}, {n[2]}")
