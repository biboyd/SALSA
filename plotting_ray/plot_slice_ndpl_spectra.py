import yt
import trident
import numpy as np
from sys import argv
import h5py
# =================== #
#
#
# =================== #

def plot_slice(ds, ray, ray_h5file, field='density'):
    """

    """

    # get beginning and end of ray
    x = ray_h5file['grid']['x']
    y = ray_h5file['grid']['y']
    z = ray_h5file['grid']['z']

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

    #Create slice along ray. keep slice pointed in z-direction
    slice = SlicePlot(ds,
                      norm_vector,
                      field,
                      north_vector = [0, 0, 1],
                      center = center)

    slice.annotate_ray(ray)
    return slice


def plot_spect(ray, out_fname):
    """

    """

    #generate spectrum defined by inputs
    spect_gen = ray.SpectrumGenerator('COS')
    spect_gen.make_spectrum(tri_ray, lines=line_list)=

    #plot spectrum
    spect_gen.plot_spectrum(out_fname)

    #check if out name is already in .png format
    if out_fname[-4:] != ".png":
        out_fname += ".png"

    return out_fname

def plot_num_density(ray_h5file, ion_name, out_fname):
    """

    """
    #get list of num density and corresponding lengths
    num_density = list(ray_h5file['grid'][ion_p_name(ion_name)+'_number_density'])
    dl_list = list(ray_h5file['grid']['dl'])

    #convert list of dl's to list of lengths from begin of ray
    num_dls = len(dl_list)
    for i in range(1, num_dls):
        dl_list[i] += dl_list[i-1]

    #make num density plots

    fig = plt.figure()
    plt.plot(dl_list, num_density)

    plt.save(out_fname)

    return fig



def compute_col_density(ray_h5file, ion_name):
    """

    """
    #get list of num density and corresponding length
    num_density = np.array(ray_h5file['grid'][ion_p_name(ion_name)+'_number_density'])
    dl_array = np.array(ray_h5file['grid']['dl'])

    #multiply num density by its dl and sum up to get column density
    col_density = np.sum( num_density*dl_array )

    return col_density
def ion_p_name(ion_name):
    """
    convert ion species name from trident style to one that
    can be used with h5 files
    """

    # 'H I' is an exception just return H
    if ion_name == 'H I':
        return 'H'

    #split up the words in ion species name
    ion_split = ion_name.split()

    #convert num from roman numeral. subtract run b/c h5
    num = tri.from_roman(ion_name[1])-1

    #combine all the names
    outname = ion_split[0] + '_p' + num
    return outname

def open_files(ds_file_name, ray_file_name):
    ds = yt.load(ds_file_name)
    ray = yt.load(ray_file_name)
    ray_h5file = h5py.File(ray_file_name)

    return ds, ray, ray_h5file

if __name__ == '__main__':
    data_set_fname = argv[1]
    ray_fname = argv[2]
    ds, ray = open_files(data_set_fname, ray_fname)
