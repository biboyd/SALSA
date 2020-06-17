from mpi4py import MPI
import numpy as np
from sys import argv
import yt
import trident
from os import makedirs, listdir
import h5py
import astropy.units  as u

from yt.data_objects.static_output import \
    Dataset

from salsa.absorber_plotter import AbsorberPlotter
from salsa.utils.filter_definitions import parse_cut_filter

def main(ds_filename, ray_dir, i_name, out_dir, frac, cut_filters, velocity_res=10):
    #init mpi
    comm = MPI.COMM_WORLD

    #setup conditions
    line_list = []#['H I', 'O VI', 'C IV']

    if comm.rank == 0:
        center, nvec, rshift, bulk_vel = (None, None, None, None)
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


    #set up plotter settings
    mp_kwargs = dict(ion_name=i_name,
    				 absorber_fields=line_list,
                     center_gal = center,
                     north_vector = nvec,
                     redshift=rshift[0],
                     frac = frac,
                     bulk_velocity=bulk_vel,
                     plot_ice=True,
                     cut_region_filters=cut_filters,
                     use_spectacle= True,
                     plot_spectacle=True,
                     velocity_res = velocity_res,
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
    create_frames(ds_filename, my_rays, impact_parameter,out_dir=out_dir, multi_plot_kwargs=mp_kwargs)
    print(f"-------------- {comm.rank} finished----------------")


def create_frames(ds_file, rays,
                  impact_parameter=None,
                  save_multi_plot=False,
                  slice_width=None,
                  slice_height=None,
                  out_dir='./',
                  multi_plot_kwargs={}, savefig_kwargs={}):
    """
    creates a movie by combining all the plots made from the ray in ray_dir

    Parameters
    -----------
        :ds_file : str or YT Dataset
            Dataset to extract absorbers from

        :rays : list
            list of full paths to rays
        slice_width : float : define width of slice in multi plot. in units kpc
        slice_height : float : define height of slice in multi plot. in units kpc
        out_dir : string/path : path to directory where frames will be saved
        multi_plot_kwargs : dict : dictionary of arguments to pass to the multi_plot
                class

    Returns:
        none
    """
    #set file names and ion name
    if isinstance(ds_file, str):
        ds = yt.load(ds_file)
    elif isinstance(ds_file, Dataset):
        ds = ds_file

    for ray_fname in rays:
        #load new multi plot and ray
        mp = AbsorberPlotter(ds, ray_filename=ray_fname, **multi_plot_kwargs)
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
            mp.fig.savefig(f"{out_dir}/plots{ray_num}.png", **savefig_kwargs)

        ray_index = int(ray_num)
        # save absorber data
        if mp.ice_df is not None:
            mp.ice_df['absorber_index']=np.empty(mp.num_ice, dtype=str)
            # add absorber index
            start = 65 # Ascii number for 'A'
            for i in range(mp.num_ice):
                mp.ice_df.loc[i,'absorber_index'] = f"{ray_num}{chr(start+i)}"

            outfile=f"{out_dir}/ray{ray_num}_ice_absorbers{mp.num_ice}.h5"
            mp.ice_df.to_hdf(outfile, mode='w')

        if mp.spectacle_df is not None:
            mp.spectacle_df['absorber_index']=np.empty(mp.num_spectacle, dtype=str)
            # add absorber index
            start = 65 # Ascii number for 'A'
            for i in range(mp.num_spectacle):
                mp.spectacle_df[i, 'absorber_index'] = f"{ray_num}{chr(start+i)}"

            outfile=f"{out_dir}/ray{ray_num}_spectacle_absorbers{mp.num_spectacle}.h5"
            mp.spectacle_df.to_hdf(outfile, mode='w')

        # close files/figures
        mp.close()


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
        cuts = argv[6] # Ex. "hot inflow cgm"
        cut_dir = "_".join(cuts.split(" "))
        out_dir +=f"/{cut_dir}"
        #retrieve properly formatted cut argument
        cut_filters = parse_cut_filter(cuts)

    elif len(argv) == 6:
        filename = argv[1]
        ray_dir = argv[2]
        ion_name = argv[3]
        out_dir= argv[4]
        frac = float(argv[5])
        cut_filters=None
        out_dir+="/no_cuts"
        print("Not applying any sort of cuts")
    else:
        raise RuntimeError("Takes 5 or 6 arguments: Dataset_fname Ray_directory Ion_name Output_directory frac (cuts)")


    #make sure out directory exists
    makedirs(out_dir, exist_ok=True)
    #check to see if should use bulk velocity
    main(filename, ray_dir, ion_name, out_dir, frac, cut_filters=cut_filters)
