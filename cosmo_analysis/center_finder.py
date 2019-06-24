import yt
import numpy as np
from sys import argv

#
#Find the center and orientation of galaxy
#return center coordinates and normal vector
#

def find_center(ds_fname):
    ds = yt.load(ds_fname)

    #let max density be center of gal
    v, center = ds.find_max('density')

    #create sphere around disk of galaxy
    sph_gal = ds.sphere(center, (50, 'kpc'))

    #compute angular momentum vector
    axis_rot = sph_gal.quantities.angular_momentum_vector()

    norm_vector = axis_rot/np.linalg.norm(axis_rot)

    return center.in_units('code_length'), norm_vector


if __name__ == '__main__':
    ds_filename = argv[1]
    c, n = find_center(ds_filename)

    print(c)
    print(n)
