from mpi4py import MPI
import numpy as np
from sys import argv, path
import yt
import trident
from os import makedirs, listdir
import h5py
import astropy.units  as u


path.insert(0, "/mnt/home/boydbre1/Repo/CGM/multi_plot_movie")
path.insert(0, "/home/bb/Repo/CGM/multi_plot_movie")
from center_finder import find_center
from multi_plot import multi_plot
def main(filename, ray_dir, i_name, out_dir, frac, cut_filter):
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
                                    cut_region_filter=cut_filter,
                                    use_spectacle= True,
                                    plot_spectacle=True,
                                    wavelength_width=20)


    #collect only ray files in ray_dir
    ray_files=[]
    impact_parameter=None
    for f in listdir(ray_dir):
        #check hdf5 file
        if (f[-3:] == ".h5"):
            full_name ="/".join((ray_dir, f))
            ray_files.append(full_name)

        #use normal_vector if saved in directory
        elif f == 'norm_vec.npy':
            full_name = "/".join((ray_dir, f))
            normal_vector = np.load(full_name)
            mp_kwargs['north_vector'] = normal_vector
        elif f == 'impact_parameter.npy':
            full_name = "/".join((ray_dir, f))
            impact_parameter = np.load(full_name)

    #split up rays betweeen proccesors
    ray_files_split = np.array_split(ray_files, comm.size)
    my_rays = ray_files_split[ comm.rank ]


    #calc the number_density limits
    num_density_range = get_num_density_range(i_name, comm, my_rays)
    mp_kwargs.update({'num_dense_min': num_density_range[0],
                      'num_dense_max': num_density_range[1]})

    #Log info for debugging
    print_rays = ""
    for r in my_rays:
        print_rays = f"{print_rays} {get_ray_num(r)}"
    print("my rank ", comm.rank, "my rays ", print_rays)

    #construct the frames
    create_frames(my_rays, impact_parameter,out_dir=out_dir, multi_plot_kwargs=mp_kwargs)
    print(f"-------------- {comm.rank} finished----------------")


def create_frames(rays,
                  impact_parameter=None,
                  save_multi_plot=False,
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


    # define info to collect
    absorber_head=np.array(['ray_num',
                            'spectacle',
                            'cloudy',
                            'start_interval',
                            'end_intervals',
                            'col density',
                            'impact_parameter',
                            'avg_velocity',
                            'x_location',
                            'y_location',
                            'z_location',
                            'radius',
                            'avg_density',
                            'avg_metallicity',
                            'avg_temperature',
                            'radial_velocity',
                            'vel_doppler'])
                            
    np.save(f"{out_dir}/absorber_info_header.npy", absorber_head)
    for ray_fname in rays:
        #load new multi plot and ray
        mp = multi_plot(ray_filename=ray_fname, **multi_plot_kwargs)
        if mp.data['l'].size == 0:
            continue

        ray_num = get_ray_num(ray_fname)
        if save_multi_plot:
            #create slice
            mp.create_slice(cmap='cividis')

            #add ray to title
            mp.slice.annotate_title(f"ray {ray_num}")

            #create and save frame
            mp.create_multi_plot(f"{out_dir}/mp{ray_num}.png")

        else:
            #just plot the 3 bottom graphs. No slice
            ax1 = mp.fig.add_subplot(311)
            ax2 = mp.fig.add_subplot(312)
            ax3 = mp.fig.add_subplot(313)

            mp.plot_num_density(ax_num_dense=ax1, ax_prop2=ax2)
            mp.plot_vel_space(ax=ax3)

            mp.fig.tight_layout()
            mp.fig.savefig(f"{out_dir}/plots{ray_num}.png")

        ray_index = int(ray_num)
        #save cloudy information
        interval_list, lcd_list = mp.get_iterative_cloud(coldens_fraction=mp.frac, min_logN=mp.cloud_min)
        abs_num=0
        for interval, lcd in zip(interval_list, lcd_list):
            if impact_parameter is None:
                impact = np.nan
            else:
                impact = impact_parameter[ray_index]
            #Create np array to store absorber data
            absorber_info = np.empty(len(absorber_head), dtype=np.float64)
            absorber_info[0] = ray_index
            absorber_info[1] = 0.
            absorber_info[2] = 1.
            absorber_info[3] = interval[0]
            absorber_info[4] = interval[1]
            absorber_info[5] = lcd
            absorber_info[6] = impact

            absorber_info[7:-1] = calc_cloudy_absorber_props(mp.ray, interval[0], interval[1])
            absorber_info[-1] = np.nan
            np.save(f"{out_dir}/ray{ray_num}_cloudy_absorber{abs_num:02d}.npy", absorber_info)
            abs_num+=1

        # save spectacle information
        if mp.spectacle_model is not None:
            lcd, del_vel, vel_doppler = calc_spectacle_absorber_props(mp.spectacle_model)
            abs_num=0
            for i in range(lcd.size):
                absorber_info = np.empty(len(absorber_head), dtype=np.float64)
                absorber_info[0] =ray_index 
                absorber_info[1] = 1.
                absorber_info[2] = 0.
                absorber_info[3] = np.nan
                absorber_info[4] = np.nan
                absorber_info[5] = lcd[i]
                absorber_info[6] = impact
                absorber_info[7] = del_vel[i]
                absorber_info[8:-1] = np.nan
                absorber_info[-1] = vel_doppler[i]

                np.save(f"{out_dir}/ray{ray_num}_spectacle_absorber{abs_num:02d}.npy", absorber_info)
                abs_num+=1


        # close files/figures
        mp.close()

def calc_spectacle_absorber_props(spec_model):
    """
    Calculate/gather data from spectacle model including absorber inforamtion
    like column density, velocity, etc.

    Parameters:
        spectacle model
    Returns:
        list of absorbers with properties
    """
    vel_array = np.linspace(-1500, 1500, 1000)*u.Unit('km/s')
    line_stats = spec_model.line_stats(vel_array)

    lcd = line_stats['col_dens']
    delta_v = line_stats['delta_v'].value
    v_doppler = line_stats['v_dop'].value


    return lcd, delta_v, v_doppler




def calc_cloudy_absorber_props(ray, start, end):
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
    props = ('velocity_los', 'x', 'y', 'z','radius', 'density', 'metallicity', 'temperature', 'radial_velocity')
    props_units=('km/s', 'kpc', 'kpc', 'kpc','kpc', 'g/cm**3', 'Zsun','K', 'km/s') 
    avg_props = []
    for prop, units in zip(props, props_units):
        #compute weighted sum of property
        avg = np.sum(dl*density*ray.data[prop][start:end])/col_dense
        avg_props.append(avg.in_units(units))

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

def get_num_density_range(ion_name, comm, my_rays, use_defaults=True):
    """
    Calculates the number density limits to be used for plotting. Either takes
    predefined limits, or calculates them by taking the limits to be:
     median/100 to median*1000

    Parameters:
        ion_name : string: name of the ion (ie H I, C IV)
        comm : MPI thing : the comm controller for the parallezation
        my_rays : list : the list of ray file paths
        use_defaults :bool : if set to True, then defaults will be used (see dict)
            instead of calculating the median value.
    Returns:
        num_density_range : numpy array: minimum and maximum limits of number density
    """
    num_density_dict = {'H I':[1e-11, 1e-5],
                       'C IV':[1e-12, 1e-6],
                       'O VI':[1e-12, 1e-6],
                       'Si III':[1e-11, 1e-5],
                       'Si II':[1e-11, 1e-5]}

    if ion_name in num_density_dict and use_defaults==True:
        num_density_range = num_density_dict[ion_name]

    else:
        #get average median value to represent scale
        num_rays = len(my_rays)
        med =0
        for i in range(num_rays):
            #load ray
            curr_ray_file = my_rays[i]
            curr_ray= yt.load(curr_ray_file)

            #get median num density
            num_density = curr_ray.data[ f"{ion_p_num(ion_name)}" ]
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

    return num_density_range
if __name__ == '__main__':
    #take in arguments
    if len(argv) == 7:
        filename = argv[1]
        ray_dir = argv[2]
        ion_name = argv[3]
        out_dir= argv[4]
        frac = float(argv[5])
        cut_filter = argv[6]

    else:
        raise RuntimeError("Takes 6 arguments: Dataset_fname Ray_directory Ion_name Output_directory frac ")

    makedirs(out_dir, exist_ok=True)
    if cut_filter == 'cgm':
        cut_filter = "((obj[('gas', 'radius')].in_units('kpc') > 10) & \
                   (obj[('gas', 'radius')].in_units('kpc') < 200)) & \
                   ((obj[('gas', 'temperature')].in_units('K') > 1.5e4) | \
                   (obj[('gas', 'density')].in_units('g/cm**3') < 2e-26))"
    else:
        cut_filter=None

    #check to see if should use bulk velocity
    main(filename, ray_dir, ion_name, out_dir, frac, cut_filter=cut_filter)
