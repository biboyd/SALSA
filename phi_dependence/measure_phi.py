import yt
import trident
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from mpl_toolkits.axes_grid1 import AxesGrid
from sys import argv
from os import listdir, makedirs
import errno

def main(ds_fname,
         center,
         north_vector,
         impact_param,
         azimuth_angle=0,
         n_rays=50,
         ray_length=300,
         ray_dir = None,
         out_dir="./"):
    """
    Creates image of slice plot and plot measuring Column Density vs polar angle.
    Rays can either be provided or constructed depending on status of ray_dir

    Parameters:
        ds_fname : string/path : path/filename of dataset to load
        center : array like object : Center of galaxy in code_length
        impact_param : float : impact parameter of rays in kpc
        azimuth_angle : float : azimuthal angle to construct rays in radians
        n_rays : int : number of rays provided or to be constructed
        ray_length : float : length of rays in kpc
        ray_dir : string/path : path to directory containing rays. If set to None
                rays will be constructed using provided information
        out_dir : string/path : directory where image will be saved and rays constructed

    Returns:
        none
    """

    ds = yt.load(ds_fname)
    slc_field = 'H_p0_number_density'

    north_vector = ds.arr(north_vector, 'dimensionless')
    center = ds.arr(center, 'code_length').in_units('kpc')
    ray_length = ds.quan(ray_length, 'kpc')
    #find axis of rotation
    axis_vector = np.array( [-north_vector[2], 0, north_vector[0]] )
    #check not zero vector
    if axis_vector[0] == 0 and axis_vector[2] == 0:
        #switch to a non zero vector
        axis_vector = np.array( [0, -north_vector[2], north_vector[1]] )

    #rotate axis_vector by azimuthal if non-zero
    if azimuth_angle != 0:
        north_rot_vector = north_vector/np.linalg.norm(north_vector)*azimuth_ang
        az_rot = Rotation.from_rotvec(north_rot_vector)
        axis_vector = az_rot.apply(axis_vector)

    #coordinates rel to center
    rs_rel, phi_array = get_coord(north_vector, impact_param, axis_vector, n_rays)
    rs_rel = ds.arr(rs_rel, 'kpc')
    ion="H I"

    #construct rays if not already made
    if ray_dir is None:
        ray_dir = construct_rays(ds, rs_rel, center, axis_vector, ray_length, n_rays, ion, out_dir)

    ray_files = []
    for file in listdir(ray_dir):
        if file[-3:] == ".h5":
            ray_files.append(file)
    #check number of rays == number of files
    if n_rays != len(ray_files):
        n_files = len(ray_files)
        raise RuntimeError(f"Number of rays specified ({n_rays}) does not equal number of files found ({n_files})")

    ray_files = sorted(ray_files)
    #get column density from rays
    col_dense = np.empty(n_rays)
    for i in range(n_rays):
        fname="/".join((ray_dir, ray_files[i]))
        ray = yt.load(f"{fname}")
        col_dense[i] = np.sum( ray.data["H_p0_number_density"] * ray.data['dl'])
        ray.close()


    slc = yt.SlicePlot(ds, normal=axis_vector,fields =slc_field,
                        center = center,
                        north_vector=north_vector,
                        width = (200, 'kpc'))
    slc.set_cmap(slc_field, cmap='magma')
    slc.set_background_color(slc_field)
    #plot markers onto slice
    for i in range(n_rays):
        slc.annotate_marker(rs_rel[i]+center, marker='.')

    #redraw slice onto figure
    fig = plt.figure(figsize=(10, 10))
    grid = AxesGrid(fig, (0.,0.,0.5,0.5),
                    nrows_ncols = (1, 1),
                    axes_pad = 0.5,
                    label_mode = "L",
                    share_all = False,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="3%",
                    cbar_pad="0%")

    #redraw slice plot onto figure
    plot = slc.plots[slc_field]
    plot.figure = fig
    plot.axes = grid[0].axes
    plot.cax = grid.cbar_axes[0]

    slc._setup_plots()

    phi_array = np.rad2deg(phi_array)
    #plot column density
    ax = fig.add_subplot(111)
    ax.plot(phi_array, col_dense)
    ax.annotate(f'Impact Parameter: {impact_param:.0f} kpc',
                xy=(0.85, 1.05),
                xycoords='axes fraction')

    ax.set_position([0, -0.25, 0.5, 0.15])
    ax.set_yscale('log')
    ax.set_title(f'{ion} Column density vs Polar Angle')
    ax.set_ylabel('Log N (N in cm^-2)')
    ax.set_xlabel('Polar Angle (degrees)')
    fig.savefig(f"{out_dir}/slice_and_angle.png", bbox_inches='tight')


def get_coord(vector, b, axis_rotation, n):
    """
    returns an array with the proper coordinates of rays rotating
    from vector to 90 degrees. also array with corresponding angles

    Parameters:
        vector : The vector in which to start rotating from
        b : Impact parameter in kpc
        axis_rotation : vector in which to rotate aboout
        n : number of points to return

    Returns:
        coordinates : array: array of points along path of rotation
        angles : array : array of the angle rotated for each point (rel to vector)
    """

    #normalize vector
    norm_vector = vector/np.linalg.norm(vector)



    #set up array to collect coordinates
    angles = np.linspace(0, np.pi/2, n)
    coordinates = np.empty((n, 3))
    del_angle = angles[1] - angles[0]

    #set up axis of rotation to be rotation vector
    rot_vector = axis_rotation/np.linalg.norm(axis_rotation) * del_angle
    rot = Rotation.from_rotvec(rot_vector)

    curr_coord = norm_vector*b
    coordinates[0] = curr_coord
    for i in range(1, n):
        curr_coord = rot.apply(curr_coord)
        coordinates[i] = curr_coord


    return coordinates, angles

def construct_rays(ds, coordinates_rel, center, axis_vector, ray_length, n_rays, ion, out_dir):
    """
    Calculates the column density for each coordinate

    Parameters:
        ds : dataset where rays are constructed
        coordinates_rel : coordinates relative to the specified center
        center : Center of galaxy/point of choosing. in kpc
        axis_vector : axis about which points were rotated
                      rays will be constructed along this direction
        ray_lenghth : the lenght of the rays in kpc

    Returns:
        col_dense : an array of column densities for each coordinate
    """

    #create ray directory and proper padding
    ray_dir = f"{out_dir}/rays"
    makedirs(ray_dir, exist_ok=True)
    pad = np.floor( np.log10(n_rays) )
    pad = int(pad) + 1

    #define offset from coord
    offset = axis_vector/np.linalg.norm(axis_vector)*ray_length/2 + center
    for i in range(n_rays):
        ray_start = coordinates_rel[i] - offset
        ray_end = coordinates_rel[i] + offset

        ray = trident.make_simple_ray(ds, ray_start, ray_end,
                                lines=[ion],
                                data_filename = f"{ray_dir}/ray{i:0{pad}d}.h5")



    return ray_dir

if __name__ == "__main__":
    ds = argv[1]

    center = [0.5, 0.5, 0.5]
    n_vec = [0, 0, 1]

    main(ds, center, n_vec, 30, ray_dir='./rays')
