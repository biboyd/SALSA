from mpi4py import MPI
import numpy as np
from sys import argv, path
import yt
import trident
from os import makedirs, listdir
import h5py

path.insert(0, "/mnt/home/boydbre1/Repo/CGM/multi_plot_movie")
path.insert(0, "/home/bb/Repo/CGM/multi_plot_movie")
from center_finder import find_center
from multi_plot import multi_plot
def main(filename, ray_dir, i_name, out_dir, frac):
    #init mpi
    comm = MPI.COMM_WORLD

    #setup conditions
    line_list = []#['H I', 'O VI', 'C IV']

    if comm.rank == 0:
        center, nvec, rshift, bulk_vel = find_center(filename)
        center = np.array( center.in_units('code_length'), dtype=np.float64)
        nvec = np.array(nvec, dtype=np.float64)
        rshift = np.array([rshift], dtype=np.float64)
        bulk_vel = np.array(bulk_vel, dtype=np.float64)
    else:
        center= np.zeros(3, dtype=np.float64)
        nvec = np.zeros(3, dtype=np.float64)
        rshift = np.zeros(1, dtype=np.float64)
        bulk_vel = np.zeros(3, dtype=np.float64)



    #broadcast information to all processes
    comm.Barrier()
    comm.Bcast([center, MPI.DOUBLE])
    comm.Bcast([nvec, MPI.DOUBLE])
    comm.Bcast([rshift, MPI.DOUBLE])

    use_bv=False
    if use_bv is True:
        comm.Bcast([bulk_vel, MPI.DOUBLE])
    else:
        bulk_vel=None
    comm.Barrier()


    #set up multiplot settings
    mp_kwargs = dict(ds_filename=filename, ion_name=i_name,
    				    absorber_fields=line_list,
                                    center_gal = center,
                                    north_vector = nvec,
                                    redshift=rshift[0],
                                    frac = frac,
                                    bulk_velocity=bulk_vel,
                                    plot_cloud=True,
                                    use_spectacle= True,
                                    plot_spectacle=True,
                                    wavelength_width=20)


    #collect only ray files in ray_dir
    ray_files=[]
    for f in listdir(ray_dir):
        #check hdf5 file
        if (f[-3:] == ".h5"):
            full_name ="/".join((ray_dir, f))
            ray_files.append(full_name)
        elif f == 'norm_vec.npy':
            full_name = "/".join((ray_dir, f))
            normal_vector = np.load(full_name)
            mp_kwargs['north_vector'] = normal_vector

    #sort the rays
    #ray_files = sorted(ray_files)

    #split up rays betweeen proccesors
    ray_files_split = np.array_split(ray_files, comm.size)
    my_rays = ray_files_split[ comm.rank ]


    #calc the number_density limits
    num_density_dict = {'H I':[1e-11, 1e-5],
                       'C IV':[1e-12, 1e-6],
                       'O VI':[1e-12, 1e-6],
                       'Si III':[1e-11, 1e-5],
                       'Si II':[1e-11, 1e-5]}

    if i_name in num_density_dict:
        num_density_range = num_density_dict[i_name]

    else:
        #get average median value to represent scale
        num_rays = len(my_rays)
        med =0
        for i in range(num_rays):
            #load ray
            curr_ray_file = my_rays[i]
            curr_ray= yt.load(curr_ray_file)

            #get median num density
            num_density = curr_ray.data[ f"{ion_p_num(i_name)}" ]
            med += np.median(num_density)

            curr_ray.close()
        med /= num_rays
        med = np.array(med)
        sum_med = np.zeros_like(med)

        #sum all medians. divide by proccesors to get average
        comm.Barrier()
        comm.Allreduce([med, MPI.DOUBLE], [sum_med, MPI.DOUBLE], op=MPI.SUM)

        avg_med = sum_med/comm.size

        #estimate min max values to number dense plot. and markers positioning
        num_density_range = np.array( [0.01*avg_med, 1000*avg_med] , dtype=np.float64)

    mp_kwargs.update({'num_dense_min': num_density_range[0],
                      'num_dense_max': num_density_range[1]})

    #create movie frames
    print_rays = ""
    for r in my_rays:
        print_rays = f"{print_rays} {get_ray_num(r)}"
    print("my rank ", comm.rank, "my rays ", print_rays)
    create_frames(my_rays, out_dir=out_dir, multi_plot_kwargs=mp_kwargs)
    print("-------------- {} finished----------------".format(comm.rank))

def create_frames(rays,
                  save_multi_plot=True,
                  slice_width=None,
                  slice_height=None,
                  out_dir='./',
                  multi_plot_kwargs={}):
    """
    creates a movie by combining all the plots made from the ray in ray_dir

    Parameters:
        rays : list : list of full paths to rays
        slice_width : float : define width of slice in multi plot. in units kpc
        slice_height : float : define height of slice in multi plot. in units kpc
        out_dir : string/path : path to directory where frames will be saved
        multi_plot_kwargs : dict : dictionary of arguments to pass to the multi_plot
                class

    Returns:
        none
    """



    absorber_head=np.array(['ray_num',
                            'start_interval',
                            'end_intervals',
                            'col density',
                            'x_location',
                            'y_location',
                            'z_location',
                            'avg metallicity',
                            'avg velocity',
                            'avg_temperature'])
    np.save(f"{out_dir}/absorber_info_header.npy", absorber_head)
    for ray_fname in rays:
        #load new multi plot and ray
        mp = multi_plot(ray_filename=ray_fname, **multi_plot_kwargs)
        mp.ray = yt.load(ray_fname)
        if save_multi_plot:
            #create slice
            mp.create_slice(cmap='cividis')

            #add ray to title
            ray_num = get_ray_num(ray_fname)
            mp.slice.annotate_title(f"ray {ray_num}")

            #create and save frame
            mp.create_multi_plot(f"{out_dir}/mp{ray_num}.png")

            #close files/figures
            mp.close()

        #save interval information
        interval_list, lcd_list = mp.get_iterative_cloud(coldens_fraction=mp.frac, min_logN=mp.cloud_min)
        abs_num=0
        for interval, lcd in zip(interval_list, lcd_list):
            #Create np array to store absorber data
            absorber_info = np.empty(len(absorber_head), dtype=np.float64)
            absorber_info[0] = ray_num
            absorber_info[1] = interval[0]
            absorber_info[2] = interval[1]
            absorber_info[3] = lcd

            absorber_info[4:] = calc_absorber_props(mp.ray, interval[0], interval[1])

            np.save(f"{out_dir}/ray{ray_num}_absorber{abs_num:02d}.npy", absorber_info)
            abs_num+=1
def calc_absorber_props(ray, start, end):
    """
    Calculate the weighted average of a list of absorber properties
    using the total gas column density as the weight.

    Parameters:
        ray : yt.ray : lightray with the fields
        start : int: beginning of interval along light ray
        end : end: end of the intervals along the light ray
    """
    dl = ray.data['dl'][start:end]
    density = ray.data[('gas', 'density')][start:end]
    col_dense = np.sum(dl*density)
    props = ('x', 'y', 'z', 'metallicity', 'velocity_los', 'temperature')
    avg_props = []
    for prop in props:
        #compute weighted sum of property
        avg = np.sum(dl*density*ray.data[prop][start:end])/col_dense
        avg_props.append(avg)

    return avg_props

def get_ray_num(file_path):
    """
    extract the ray's number from it's file name

    Parameters:
        file_path : string: full path to ray file
    Returns:
        num :string : the number corresponding to the ray
    """
    filename = file_path.split('/')[-1]
    num = filename[3:-3]

    return num
def ion_p_num(ion_name):
    """
    convert ion species name from trident style to one that
    can be used with h5 files
    """
    #split up the words in ion species name
    ion_split = ion_name.split()
    #convert num from roman numeral. subtract run b/c h5
    num = trident.from_roman(ion_split[1])-1

    #combine all the names
    outname = f"{ion_split[0]}_p{num}_number_density"
    return outname

if __name__ == '__main__':
    #take in arguments
    if len(argv) == 6:
        filename = argv[1]
        ray_dir = argv[2]
        ion_name = argv[3]
        out_dir= argv[4]
        frac = float(argv[5])

    else:
        raise RuntimeError("Takes 5 arguments: Dataset_fname Ray_directory Ion_name Output_directory frac ")

    #check to see if should use bulk velocity
    main(filename, ray_dir, ion_name, out_dir, frac)
