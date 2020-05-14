import yt
import trident
import numpy as np
from mpi4py import MPI
from sys import argv
from os import makedirs
from scipy.spatial.transform import Rotation

from CGM.general_utils.construct_rays import construct_rays
from CGM.general_utils.center_finder import find_center
from CGM.general_utils.filter_definitions import radius_function

def movie_sightlines(dataset,
                    length=200,
                    n_rays=100,
                    norm_vector=[0, 0, 1],
                    az_angle=0,
                    off_angle=45,
                    max_impact_param=200,
                    center=None):

    """
    Parameters:
        dataset: enzo dataset on which to construct the rays
        length: the length of the rays in kpc
        n_rays: the number of rays to construct
        norm_vector: normal vector of galaxy. (perpendiuclar to disk)
        az_angle: The azimuthal angle around the direction. in degrees
        off_angle: Polar angle to normal_vector. in degrees.
        max_impact_param: sets range to construct rays. in kpc from center of galaxy
        center : array/list. Coordinates to center ray creation. Defaults to dataset/s center. in kpc

    Returns:
        start_points : numpy array : the starting coordinates for each lightray
        end_points : numpy array : the ending coordinates for each lightray
        norm_vecotr : numpy array : normal vector that will fix galaxy orientation
                when plotting using multi_plot class.

    """

    ds = yt.load(dataset)
    #convert lengths to code_length
    length = ds.quan(length, 'kpc').in_units('code_length')
    max_impact_param = ds.quan(max_impact_param, 'kpc').in_units('code_length')
    norm_vector = norm_vector/np.linalg.norm(norm_vector)

    #create the unit vector for rays
    ray_unit = np.array( [norm_vector[2], 0, -norm_vector[0]] )

    #check not zero vector
    if ray_unit[0] == 0 and ray_unit[2] == 0:
        #switch to a non zero vector
        ray_unit = np.array( [0, norm_vector[2], -norm_vector[1]] )

    ray_unit = ray_unit/np.linalg.norm(ray_unit)
    #rotate ray unit
    if az_angle != 0:
        angle = np.deg2rad(az_angle)
        rot_vector = norm_vector*angle
        rot = Rotation.from_rotvec(rot_vector)
        ray_unit = rot.apply(ray_unit)

    if off_angle != 0:
        #create rot vector perpendicular to normal and ray unit vec
        rot_vector = np.cross(ray_unit, norm_vector)

        #convert degree to rad and scale rot_vector
        off_rad = np.deg2rad(off_angle)
        rot_vector *= off_rad/np.linalg.norm(rot_vector)

        #apply rotation to ray unit and norm vectors
        rot= Rotation.from_rotvec(rot_vector)
        ray_unit = rot.apply(ray_unit)
        norm_vector = rot.apply(norm_vector)

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
    start_points = ray_centers - offset
    end_points = ray_centers + offset

    return start_points, end_points, norm_vector

def movie_rays( dataset,
                    line_list,
                    length=200,
                    n_rays=100,
                    norm_vector=[0, 0, 1],
                    az_angle=0,
                    off_angle=45,
                    bulk_vel = None,
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
        az_angle: The azimuthal angle around the direction. in degrees
        off_angle: Polar angle to normal_vector. in degrees.
        max_impact_param: sets range to construct rays. in kpc from center of galaxy
        center : array/list. Coordinates to center ray creation. Defaults to dataset/s center. in kpc
        parallel : boolean. If parallelization is being used set to true.
        out_dir: directory in which to save the rays

    Returns:
        none
    """
    #create directory if doesn't exist
    makedirs(out_dir, exist_ok=True)

    other_fields=['density', 'metallicity', 'temperature', ('gas', 'radius'), 'radial_velocity']
    ds = yt.load(dataset)
    #add ion fields to dataset if not already there
    trident.add_ion_fields(ds, ions=line_list, ftype='gas')

    # add radius field to dataset
    ds.add_field(('gas', 'radius'),
             function=radius_function,
             units="code_length",
             take_log=False,
             validators=[yt.fields.api.ValidateParameter(['center'])])

    #collect sightlines
    start_points, end_points, norm_vector = movie_sightlines(dataset,
                                                length=length,
                                                n_rays=n_rays,
                                                norm_vector=norm_vector,
                                                az_angle=az_angle,
                                                off_angle=off_angle,
                                                max_impact_param=max_impact_param,
                                                center=center)

    #save normal vector for future plotting purposes
    np.save(f"{out_dir}/norm_vec.npy", norm_vector)

    #construct lightrays
    construct_rays(ds, start_points, end_points,
                   center=center, bulk_velocity=bulk_vel,
                   line_list=line_list,
                   other_fields=other_fields,
                   out_dir=out_dir,
                   parallel=parallel)

if __name__ == '__main__':
    #setup conditions
    line_list = ['H I','H II','Si II', 'Si III', 'C IV', 'O VI', 'Ne VIII', 'Mg X']
    if len(argv) == 5:
        filename = argv[1]
        num_rays=int(argv[2])
        ray_length=int(argv[3])
        out_dir = argv[4]
    else:
        raise RuntimeError("Takes in 4 Arguments. Dataset_filename num_rays ray_lenght out_directory")

    center, n_vec, rshift, bv = find_center(filename)
    #divide rays evenly
    movie_rays(filename, line_list,
                    n_rays=num_rays,
                    norm_vector = n_vec,
                    bulk_vel = bv,
                    length=ray_length,
                    max_impact_param=ray_length/2,
                    center = center,
                    out_dir=out_dir,
                    parallel = True)
