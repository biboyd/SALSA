import yt
import trident
import numpy as np
from mpi4py import MPI
from sys import argv
from os import makedirs
import errno


def construct_rays( dataset,
                    line_list,
                    length=200,
                    n_rays=100,
                    direction='z',
                    angle=0,
                    dist_from_center=200,
                    center=None,
                    parallel=False,
                    out_dir='./rays',
                    DEBUG=False):
    """
    Constructs a number of light rays to "scan" a galactic data set using trident

    Parameters:
        dataset: enzo dataset on which to construct the rays
        line_list: list of ions to include in rays
        length: the length of the rays in kpc
        n_rays: the number of rays to construct
        direction: the coordinate direction in which to "scan" the galaxy
        angle: The azimuthal angle around the direction. in degrees
        dist_from_center: range to construct rays. in kpc from center of galaxy
        center : array/list. Coordinates to center ray creation. Defaults to dataset/s center. in kpc
        parallel : boolean. If parallelization is being used set to true.
        out_dir: directory in which to save the rays

    Returns:
        none
    """
    if DEBUG:
        print(trident.__version__)
    ds = yt.load(dataset)
    #add ion fields to dataset if not already there
    trident.add_ion_fields(ds, ions=line_list, ftype='gas')

    #create directory if doesn't exist
    try:
        makedirs(out_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    #start MPI in case parallelization is being used
    if parallel:
        comm = MPI.COMM_WORLD

    #convert lengths to code_length and angle to radians
    length = ds.quan(length, 'kpc').in_units('code_length')
    dist_from_center = ds.quan(dist_from_center, 'kpc').in_units('code_length')
    angle = np.deg2rad(angle)

    if center is None:
        center = ds.domain_center.in_units('code_length')
    else:
        center = ds.arr(center, 'kpc').in_units('code_length')

    #get right indexing for direction
    if (direction == 'x'):
        dir_index=0
        coord1_index=1
        coord2_index=2
    elif (direction == 'y'):
        dir_index=1
        coord1_index=2
        coord2_index=0
    elif (direction == 'z'):
        dir_index=2
        coord1_index=0
        coord2_index=1
    else:
        raise RuntimeError("direction must be 'x', 'y', or 'z'")

    #calculate the changing coordinate
    begin = center[dir_index] + dist_from_center
    end = center[dir_index] - dist_from_center
    direct_coord = np.linspace(begin, end, n_rays)

    #compute ray length in the coordinate directions
    len_coord1 = length/2* np.cos(angle)
    len_coord2 = length/2* np.sin(angle)

    #calc beginning and ending of ray for constant coordinates
    coord1_begin = center[coord1_index] - len_coord1
    coord1_end = center[coord1_index] + len_coord1

    coord2_begin = center[coord2_index] - len_coord2
    coord2_end = center[coord2_index] + len_coord2

    #empty ray coordinates to be filled
    ray_begin = ds.arr(np.zeros(3), "code_length")
    ray_end = ds.arr(np.zeros(3), "code_length")

    #set padding for filenames
    pad = np.floor( np.log10(n_rays) )
    pad = int(pad) + 1


    my_ray_nums = np.arange(n_rays)
    if parallel:
        #split ray numbers then take a portion based on rank
        split_ray_nums = np.array_split(my_ray_nums, comm.size)
        my_ray_nums = split_ray_nums[ comm.rank ]
    for i in my_ray_nums:
        #set beginning ray
        ray_begin[dir_index] = direct_coord[i]
        ray_begin[coord1_index] = coord1_begin
        ray_begin[coord2_index] = coord2_begin

        #set ending ray
        ray_end[dir_index] = direct_coord[i]
        ray_end[coord1_index] = coord1_end
        ray_end[coord2_index] = coord2_end

        #construct ray
        trident.make_simple_ray(ds,
                                ray_begin,
                                ray_end,
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
        n_rays=int(argv[2])
        ray_length=int(argv[3])
        out_dir = argv[4]
    else:
        raise RuntimeError("Takes in 4 Arguments. Dataset_filename num_rays ray_lenght out_directory")


    #divide rays evenly
    construct_rays(filename, line_list,
                    n_rays=n_rays,
                    length=ray_length,
                    dist_from_center=ray_length/2,
                    out_dir=out_dir,
                    parallel = True)

