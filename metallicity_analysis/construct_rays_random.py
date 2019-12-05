import yt
import trident
import numpy as np

from mpi4py import MPI
from sys import argv
from os import makedirs
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt
from sys import path
path.insert(0, "/mnt/home/boydbre1/Repo/CGM/multi_plot_movie")
path.insert(0, "/home/bb/Repo/CGM/multi_plot_movie")
from center_finder import find_center
#sample impact param

def random_sightlines(dsname, center, num_sightlines, max_impact_param, min_impact_param=0, length=200, seed=None):
    """
    randomly sample impact parameter to get random sightlines from a give galaxy

    Parameters:
        center : arr: coordinates of the center of the galaxy
        num_sightlines : int : number of sightlines to return
        max_impact_param : float : maximum impact param to sample from in kpc
        min_impact_param : float : minimum impact param to sample from in kpc
        length : float : length of the sightline in kpc

    Returns :
        start_points : array : 2d array of the startpoints for
                    each sightline
        end_points : array : 2d array of the endpoints for
                    each sightline
    """

    np.random.seed(seed)
    ds = yt.load(dsname)
    length = ds.quan(length, 'kpc').in_units('code_length')
    length = length.value
    min_impact_param = ds.quan(min_impact_param, 'kpc').in_units('code_length')
    max_impact_param = ds.quan(max_impact_param, 'kpc').in_units('code_length')


    #take sqrt so that impact param is uniform in projected area space
    impact_param = np.sqrt(np.random.uniform(min_impact_param.value**2, max_impact_param.value**2, num_sightlines))
    #theta represents polar angle. phi represents azimuthal
    theta = np.random.uniform(0, np.pi, num_sightlines)
    phi = np.random.uniform(0, 2*np.pi, num_sightlines)

    rad_vec= np.empty((num_sightlines, 3))
    rad_vec[:, 0] = impact_param*np.cos(phi)*np.sin(theta)
    rad_vec[:, 1] = impact_param*np.sin(phi)*np.sin(theta)
    rad_vec[:, 2] = impact_param*np.cos(theta)

    #define vector along sightline (perpendicular to radial vector)
    perp_vec = np.empty_like(rad_vec)
    perp_vec[:, 0] = rad_vec[:, 1]
    perp_vec[:, 1] = -1*rad_vec[:, 0]
    perp_vec[:, 2] = 0.

    #randomly rotate perp_vec around rad vec
    alpha=np.random.uniform(0., 2*np.pi, num_sightlines)
    for i in range(num_sightlines):
        #normalize perpendicular vector
        perp_vec[i, :] =  perp_vec[i, :]/np.sqrt(perp_vec[i, 0]**2 + perp_vec[i, 1]**2)
        #rotate around radial vector by random amount
        rot_vec = alpha[i] *rad_vec[i, :]/np.linalg.norm(rad_vec[i, :])
        rot = Rotation.from_rotvec(rot_vec)
        perp_vec[i, :] = rot.apply(perp_vec[i, :])


    #shift to be centered at galaxy
    sightline_centers = rad_vec +np.array(center)
    #find ending and start points for each sightline
    end_point = sightline_centers + length/2 *perp_vec
    start_point = sightline_centers - length/2 *perp_vec

    return start_point, end_point

def random_rays(dsname, center,
                n_rays, max_impact,
                min_impact_param=0.,
                length=200,
                line_list=['H I', 'C IV', 'O VI'],
                other_fields=['density', 'metallicity'],
                out_dir='./',
                parallel=True,
                seed=None):
    """

    """
    #stuff
    if parallel:
        comm = MPI.COMM_WORLD

    ds = yt.load(dsname)
    start_points, end_points = random_sightlines(dsname, center,
                                                 n_rays,
                                                 max_impact,
                                                 min_impact_param=min_impact_param,
                                                 length=length,
                                                 seed=seed)

    #set padding for filenames
    pad = np.floor( np.log10(n_rays) )
    pad = int(pad) + 1

    my_ray_nums = np.arange(n_rays)
    if parallel:
        #split ray numbers then take a portion based on rank
        split_ray_nums = np.array_split(my_ray_nums, comm.size)
        my_ray_nums = split_ray_nums[ comm.rank]
    for i in my_ray_nums:
        #construct ray
        trident.make_simple_ray(ds,
                                start_points[i],
                                end_points[i],
                                lines=line_list, fields=other_fields,
                                data_filename= f"{out_dir}/ray{i:0{pad}d}.h5")
if __name__ == '__main__':
    #setup conditions
    line_list = ['H I','H II','Si II', 'Si III', 'C IV', 'O VI', 'Ne VIII', 'Mg X']
    if len(argv) == 7:
        filename = argv[1]
        num_rays=int(argv[2])
        ray_length=int(argv[3])
        min_impact= int(argv[4])
        max_impact=int(argv[5])
        out_dir = argv[6]
    else:
        raise RuntimeError("Takes in 5 Arguments. Dataset_filename num_rays ray_lenght max_impact_param out_directory")

    center, n_vec, rshift, bv = find_center(filename)
    random_rays(filename, center, num_rays, max_impact, min_impact_param=min_impact, line_list=line_list, length=ray_length, out_dir=out_dir)
