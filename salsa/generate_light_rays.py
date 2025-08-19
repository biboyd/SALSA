import os
import yt
import trident
import numpy as np
import h5py

from yt.data_objects.static_output import \
    Dataset

from mpi4py import MPI

from scipy.spatial.transform import Rotation

from unyt.array import unyt_quantity
from astropy.units.quantity import Quantity

def random_sightlines(ds_file,
                      center,
                      num_sightlines, 
                      max_impact_param,
                      min_impact_param=0,
                      length=200):
    """
    randomly sample impact parameter to get random sightlines from a given galaxy center

    Parameters
    ----------
    ds_file : str or YT dataset
        path to dataset to loar or already loaded dataset

    center : array like
        coordinates of the center of the galaxy in units code_length

    num_sightlines : int
        number of sightlines to return

    max_impact_param : float
        maximum impact param to sample from. If no units are attached
        (either unyt or astropy) it is assumed to be in kpc

    min_impact_param : float, optional
        minimum impact param to sample from. If no units are attached
        (either unyt or astropy) it is assumed to be in kpc
        Default: 0.

    length : float, optional
         length of the sightline. f no units are attached
        (either unyt or astropy) it is assumed to be in kpc
         Default: 200

    Returns
    --------
    start_points : array
        2d array of the startpoints for each sightline in code_length

    end_points : array
        2d array of the endpoints for each sightline in code_length

    impact_param : array
        array of impact parameters for each ray created in code_length
    """

    #set file names and ion name
    if isinstance(ds_file, str):
        ds = yt.load(ds_file)
    elif isinstance(ds_file, Dataset):
        ds = ds_file

    # check for units (unyt or astropy) and convert to unyt/yt code_length
    if isinstance(min_impact_param, unyt_quantity):
        min_impact_param = ds.quan(min_impact_param).in_units('code_length')
    elif isinstance(min_impact_param, Quantity):
        min_impact_param = ds.quan(min_impact_param, 
                                   str(min_impact_param.unit)).in_units('code_length')
    else:
        min_impact_param = ds.quan(min_impact_param, 'kpc').in_units('code_length')

    if isinstance(max_impact_param, unyt_quantity):
        max_impact_param = ds.quan(max_impact_param).in_units('code_length')
    elif isinstance(max_impact_param, Quantity):
        max_impact_param = ds.quan(max_impact_param, 
                                   str(max_impact_param.unit)).in_units('code_length')
    else:
        max_impact_param = ds.quan(max_impact_param, 'kpc').in_units('code_length')

    if isinstance(length, unyt_quantity):
        length = ds.quan(length).in_units('code_length')
    elif isinstance(length, Quantity):
        length = ds.quan(length, 
                         str(length.unit)).in_units('code_length')
    else:
        length = ds.quan(length, 'kpc').in_units('code_length')
    length = length.value

    #randomly select angle and distance from center of gal
    #take sqrt so that impact param is uniform in projected area space
    impact_param = np.sqrt(np.random.uniform(min_impact_param**2,
                                             max_impact_param**2,
                                             num_sightlines))

    #theta represents polar angle. phi represents azimuthal
    theta = np.random.uniform(0, np.pi, num_sightlines)
    phi = np.random.uniform(0, 2*np.pi, num_sightlines)

    #construct vector from gal_center to sightline midpoint
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

    return start_point, end_point, impact_param

def construct_rays(ds_file,
        start_points,
        end_points,
        impact_params,
        fld_params=None,
        line_list=None,
        other_fields=None,
        ftype='gas',
        ray_directory='./'):
    """
    Construct rays given a set of starting points and end points.

    Parameters
    ----------
    ds_file : str or YT dataset
        path to dataset to be used to create rays

    start_points : numpy array
        1d array of starting points for each ray (code_length)

    end_points : numpy array
        1d array of end points for each ray (code_length)

    impact_params : numpy array
        1d array of impact parameters for each ray (code_length).
        Used to save the impact parameter to the ray dataset.

    fld_params: dict, optional
        Dictionary of parameters that will be passed to the lightrays. (ie
        `center`, `bulk_velocity`).
        Default: None

    line_list : list
        list of ions to add to light rays. None defaults to
        H I, C IV, and O VI

    other_fields : list
        other yt fields to add to light rays. None defaults
        to density, metallicity, and temperature

    ftype : str
        The field to be passed to trident that ion fields will be added to, i.e.
        ('gas', 'H_p0_number_density'). 'gas' should work for most grid-based
        simulations. For particle-based simulations this will not work and needs
        to be changed. 'PartType0' often works though it varies.
        See trident.make_simple_ray() for more information

    ray_directory : str/path
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

    if comm.rank == 0:
        os.makedirs(ray_directory, exist_ok=True)
    comm.Barrier()
    
    for i in my_ray_nums:
        #construct ray
        ray_filename = f"{ray_directory}/ray{i:0{pad}d}.h5"
        yt_ray = trident.make_simple_ray(ds_file,
                                         start_points[i],
                                         end_points[i],
                                         lines=line_list,
                                         fields=other_fields,
                                         ftype=ftype,
                                         field_parameters=fld_params,
                                         data_filename=ray_filename)
        yt_ray.close()  # we no longer need the yt object in memory
                        # even though this function doesn't actually do anything...

        with h5py.File(ray_filename, "a") as f:
            f.attrs["light_ray_solution_impact_parameter"] = [impact_params[i]]
            f.attrs["light_ray_solution_impact_parameter_units"] = "code_length"

    comm.Barrier()

def generate_lrays(ds, center,
                n_rays, max_impact_param,
                min_impact_param=0.,
                length=200,
                fld_params={},
                ion_list=['H I', 'C IV', 'O VI'],
                fields=None,
                ftype='gas',
                ray_directory='./'):
    """
    Generate a sample of trident lightrays that randomly, uniformly cover
    impact parameter.

    Parameters
    ----------
    ds_file : str or YT Dataset
        path to dataset or already loaded, YT dataset

    center : arr
        coordinates of the center of the galaxy

    n_rays : int
        number of light rays to construct

    max_impact_param : float
        maximum impact param to sample from in kpc

    min_impact_param : float
        minimum impact param to sample from in kpc

    length : float
        length of the sightline in kpc

    ion_list : list
        ions to add to lightray

    fields : list
        fields to add to lightray

    ftype : str
        The field to be passed to trident that ion fields will be added to, i.e.
        ('gas', 'H_p0_number_density'). 'gas' should work for most grid-based
        simulations. For particle-based simulations this will not work and needs
        to be changed. 'PartType0' often works though it varies.
        See trident.add_ion_fields() for more information

    ray_directory : string
        path to where ray files will be written

    """

    comm = MPI.COMM_WORLD
    #collect sightlines
    if comm.rank == 0:
        start_pnts, end_pnts, imp_param = random_sightlines(ds, center,
                                                 n_rays,
                                                 max_impact_param,
                                                 min_impact_param=min_impact_param,
                                                 length=length)
        
        # os.makedirs(ray_directory, exist_ok=True)
        # imp_param = ds.arr(imp_param, 'code_length').in_units('kpc')
        # np.save(f"{ray_directory}/impact_parameter.npy", imp_param)

    else:
        start_pnts= np.empty((n_rays, 3), dtype=np.float64)
        end_pnts= np.empty((n_rays, 3), dtype=np.float64)
        imp_param= np.empty(n_rays, dtype=np.float64)

    # share sightline points
    comm.Barrier()
    comm.Bcast([start_pnts, MPI.DOUBLE])
    comm.Bcast([end_pnts, MPI.DOUBLE])
    comm.Bcast([imp_param, MPI.DOUBLE])

    #add center to field parameters
    fld_params['center']=center

    #add density field. Needed in absorber calculations
    if fields is None:
        construct_fields=['density']
    else:
        construct_fields = fields.copy()
        if 'density' not in construct_fields:
            construct_fields.append('density')

        # remove x/y/z b/c will cause problems with trident
        for coord in ['x', 'y', 'z']:
            if coord in construct_fields:
                construct_fields.remove(coord)

    #add ion fields to dataset if not already there
    trident.add_ion_fields(ds, ions=ion_list, ftype=ftype)

    construct_rays(ds, start_pnts, end_pnts, imp_param,
                   fld_params=fld_params,
                   line_list=ion_list,
                   other_fields=construct_fields,
                   ftype=ftype,
                   ray_directory=ray_directory)
