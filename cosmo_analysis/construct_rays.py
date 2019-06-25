import yt
import trident
import numpy as np
from center_finder import find_center
from mpi4py import MPI
from sys import argv
from os import makedirs
import errno
from scipy.spatial.transform import Rotation


def construct_rays( dataset,
                    line_list,
                    length=200,
                    n_rays=100,
                    norm_vector=[0, 0, 1],
                    angle=0,
                    max_impact_param=200,
                    center=None,
                    parallel=False,
                    out_dir='./rays'):
    """
    Constructs a number of light rays to "scan" a galactic data set using trident

    Parameters:
        dataset: enzo dataset on which to construct the rays
        line_list: list of ions to include in rays
        length: the length of the rays in kpc
        n_rays: the number of rays to construct
        norm_vector: normal vector of galaxy. (perpendiuclar to disk)
        angle: The azimuthal angle around the direction. in degrees
        max_impact_param: sets range to construct rays. in kpc from center of galaxy
        center : array/list. Coordinates to center ray creation. Defaults to dataset/s center. in kpc
        parallel : boolean. If parallelization is being used set to true.
        out_dir: directory in which to save the rays

    Returns:
        none
    """

    ds = yt.load(dataset)
    #add ion fields to dataset if not already there
    trident.add_ion_fields(ds, ions=line_list, ftype='gas')

    #create directory if doesn't exist
    makedirs(out_dir, exist_ok=True)


    #start MPI in case parallelization is being used
    if parallel:
        comm = MPI.COMM_WORLD

    #convert lengths to code_length and angle to radians
    length = ds.quan(length, 'kpc').in_units('code_length')
    max_impact_param = ds.quan(max_impact_param, 'kpc').in_units('code_length')
    norm_vector = norm_vector/np.linalg.norm(norm_vector)

    #create the unit vector for rays
    ray_unit = np.array( [norm_vector[2], 0, -norm_vector[0]] )
    ray_unit = ray_unit/np.linalg.norm(ray_unit)

    #check not zero vector
    if ray_unit[0] == 0 and ray_unit[2] == 0:
        #switch to a non zero vector
        ray_unit = np.array( [0, norm_vector[2], -norm_vector[1]] )

    #rotate ray unit
    if angle != 0:
        angle = np.deg2rad(angle)
        rot_vector = norm_vector*angle
        rot = Rotation.from_rotvec(rot_vector)
        ray_unit = rot.apply(ray_unit)

    if center is None:
        center = ds.domain_center.in_units('code_length')
    else:
        center = ds.arr(center, 'code_length')

    #find the beginning and ending centers of all rays
    start_ray_cent = center + max_impact_param*norm_vector
    end_ray_cent = center - max_impact_param*norm_vector
    ray_centers = np.linspace(start_ray_cent, end_ray_cent, n_rays)

    #add offset to find end and beginning to all rays
    offset = length/2 * ray_unit
    ray_begins = ray_centers - offset
    ray_ends = ray_centers + offset

    #set padding for filenames
    pad = np.floor( np.log10(n_rays) )
    pad = int(pad) + 1

    my_ray_nums = np.arange(n_rays)
    if parallel:
        #split ray numbers then take a portion based on rank
        split_ray_nums = np.array_split(my_ray_nums, comm.size)
        my_ray_nums = split_ray_nums[ comm.rank ]

    for i in my_ray_nums:
        #construct ray
        trident.make_simple_ray(ds,
                                ray_begins[i],
                                ray_ends[i],
                                lines=line_list,
                                data_filename= f"{out_dir}/ray{i:0{pad}d}.h5")
    if parallel:
        comm.Barrier()
        if comm.rank == 0:
             print("-----all finished------")

#now actual test
if __name__ == '__main__':
    #setup conditions
    line_list = ['H I', 'O VI', 'C IV']
    if len(argv) == 5:
        filename = argv[1]
        num_rays=int(argv[2])
        ray_length=int(argv[3])
        out_dir = argv[4]
    else:
        raise RuntimeError("Takes in 4 Arguments. Dataset_filename num_rays ray_lenght out_directory")

    center, n_vec = find_center(filename)
    #divide rays evenly
    construct_rays(filename, line_list,
                    n_rays=num_rays,
                    norm_vector = n_vec,
                    length=ray_length,
                    max_impact_param=ray_length/2,
                    center = center,
                    out_dir=out_dir,
                    parallel = True)
