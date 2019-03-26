import yt
import trident
import numpy as np
from sys import argv
import h5py
# =================== #
#
#
#

def plot_slice(ds_file_name, ray_file_name, field='density'):
    """

    """
    ds = yt.load(ds_file_name)
    ray = yt.load(ray_file_name)

    # get beginning and end of ray
    f = h5py.File(ray_file_name)
    x = f['grid']['x']
    y = f['grid']['y']
    z = f['grid']['z']

    ray_begin = np.array([ x[0], y[0], z[0] ])
    ray_end = np.array([ x[-1], y[-1], z[-1] ])

    #construct vector pointing in ray's direction
    ray_vec = ray_end - ray_begin

    #construct vec orthogonal to ray
    norm_vector = [ray[1], -1*ray[0], 0]

    #get center of ray keep z at center
    center = ds.domain_center

    ray_center = (ray_begin + ray_end)/2
    ray_center[2] = center[2]

    #Create slice along ray. keep image pointed in z-dir
    slice = SlicePlot(ds,
                      norm_vector,
                      field,
                      north_vector = [0, 0, 1],
                      center = center)

    slice.annotate_ray(ray)
    return slice


if __name__ == '__main__':
    data_set_fname = argv[1]
    ray_fname = argv[2]
