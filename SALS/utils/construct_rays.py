import yt
import trident
import numpy as np

from mpi4py import MPI

def construct_rays(ds_file,
        start_points,
        end_points,
        fld_params=None,
        line_list=None,
        other_fields=None,
        out_dir='./'):
    """
    Construct rays using trident.

    Parameters
    ----------
        :ds_file : str or YT dataset
            path to dataset to be used to create rays

        :start_points : numpy array
            1d array of starting points for each ray (code_length)

        :end_points : numpy array
            1d array of end points for each ray (code_length)

        :fld_params: dict, optional
            Dictionary of parameters that will be passed to the lightrays. (ie
            `center`, `bulk_velocity`).
            Default: None

        :line_list : list
            list of ions to add to light rays. None defaults to
            H I, C IV, and O VI

        :other_fields : list
            other yt fields to add to light rays. None defaults
            to density, metallicity, and temperature

        :out_dir : str/path
            where to save all of the lightrays
    """
    comm = MPI.COMM_WORLD

    #set defaults
    if line_list is None:
        line_list=['H I', 'C IV', 'O VI']
    if other_fields is None:
        other_fields=['density', 'metallicity', 'temperature']

    n_rays = start_points.shape[0]
    #set padding for filenames
    pad = np.floor( np.log10(n_rays) )
    pad = int(pad) + 1

    # distribute rays to proccesors
    my_ray_nums = np.arange(n_rays)

    #split ray numbers then take a portion based on rank
    split_ray_nums = np.array_split(my_ray_nums, comm.size)
    my_ray_nums = split_ray_nums[ comm.rank]

    for i in my_ray_nums:
        #construct ray
        ray_filename = f"{out_dir}/ray{i:0{pad}d}.h5"
        trident.make_simple_ray(ds_file,
                                start_points[i],
                                end_points[i],
                                lines=line_list,
                                fields=other_fields,
                                field_parameters=fld_params,
                                data_filename=ray_filename)

    comm.Barrier()
