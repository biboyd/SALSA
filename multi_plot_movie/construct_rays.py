import yt
import trident
import numpy as np
from center_finder import find_center
from mpi4py import MPI
from sys import argv
from os import makedirs
from scipy.spatial.transform import Rotation


def construct_rays( dataset,
                    line_list,
                    length=200,
                    n_rays=100,
                    norm_vector=[0, 0, 1],
                    angle=0,
                    off_angle=45,
                    bulk_vel = 0,
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
        offset: Polar angle to normal_vector. in degrees.
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

    # add radius field to dataset
    ds.add_field(('gas', 'radius'),
             function=_radius,
             units="code_length",
             take_log=False,
             validators=[yt.fields.api.ValidateParameter(['center'])])

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

    #check not zero vector
    if ray_unit[0] == 0 and ray_unit[2] == 0:
        #switch to a non zero vector
        ray_unit = np.array( [0, norm_vector[2], -norm_vector[1]] )

    ray_unit = ray_unit/np.linalg.norm(ray_unit)
    #rotate ray unit
    if angle != 0:
        angle = np.deg2rad(angle)
        rot_vector = norm_vector*angle
        rot = Rotation.from_rotvec(rot_vector)
        ray_unit = rot.apply(ray_unit)

    if off_angle != 0:
        #create rot vector perpendicular to normal and ray unit vec
        rot_vector = np.cross(ray_unit, norm_vector)

        #convert degree to rad and scale rot_vector
        offset_rad = np.deg2rad(off_angle)
        rot_vector *= offset_rad/np.linalg.norm(rot_vector)

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
    ray_begins = ray_centers - offset
    ray_ends = ray_centers + offset

    #set padding for filenames
    pad = np.floor( np.log10(n_rays) )
    pad = int(pad) + 1

    my_ray_nums = np.arange(n_rays)
    if parallel:
        #split ray numbers then take a portion based on rank
        split_ray_nums = np.array_split(my_ray_nums, comm.size)
        my_ray_nums = split_ray_nums[ comm.rank]

    #define center of galaxy for lrays
    if center is not None:
        fld_param = {'center':center}
    else:
        fld_param=None
    for i in my_ray_nums:
        #construct ray
        trident.make_simple_ray(ds,
                                ray_begins[i],
                                ray_ends[i],
                                lines=line_list,
                                fields= ['density', 'metallicity', 'temperature', ('gas', 'radius')],
                                field_parameters=fld_param,
                                data_filename= f"{out_dir}/ray{i:0{pad}d}.h5")
    if parallel:
        comm.Barrier()
        if comm.rank == 0:
             print("-----all finished------")

    return norm_vector

#function to create field in yt
def _radius(field, data):
    if data.has_field_parameter("center"):
        c = data.get_field_parameter("center")
    else:
        c = data.ds.domain_center

    x = data[('gas', 'x')] - c[0]
    y = data[('gas', 'y')] - c[1]
    z = data[('gas', 'z')] - c[2]
    return np.sqrt(x*x + y*y + z*z)


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
    new_n_vec = construct_rays(filename, line_list,
                    n_rays=num_rays,
                    norm_vector = n_vec,
                    bulk_vel = bv,
                    length=ray_length,
                    max_impact_param=ray_length/2,
                    center = center,
                    out_dir=out_dir,
                    parallel = True)
    # save vector to fix orientation for future plots 
    np.save(f"{out_dir}/norm_vec.npy", new_n_vec)
