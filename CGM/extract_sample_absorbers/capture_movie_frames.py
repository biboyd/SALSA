from mpi4py import MPI
import numpy as np
from sys import argv
import yt
import trident
from CGM.absorber_extraction_class.absorber_plotter import absorber_plotter
from CGM.general_utils.center_finder import find_center
from CGM.general_utils.filter_definitions import ion_p_num
from os import makedirs, listdir
import h5py

def main(filename, ray_dir, i_name, out_dir, use_bv, frac, cut_list=None):
    #init mpi
    comm = MPI.COMM_WORLD

    #find setup conditions
    center, nvec, rshift, bulk_vel = find_center(filename)
    center = center.in_units('code_length')

    #set up plotter settings
    mp_kwargs = dict(ds_filename=filename, ion_name=i_name,
                     cut_region_filter=cut_list,
                     center_gal = center,
                     north_vector = nvec,
                     redshift=rshift,
                     frac = frac,
                     bulk_velocity=None,
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
                  slice_width=None,
                  slice_height=None,
                  out_dir='./',
                  save_data=False,
                  multi_plot_kwargs={}):
    """
    creates a movie by combining all the plots made from the ray in ray_dir

    Parameters:
        rays : list : list of full paths to rays
        slice_width : float : define width of slice in multi plot. in units kpc
        slice_height : float : define height of slice in multi plot. in units kpc
        out_dir : string/path : path to directory where frames will be saved
        save_data : bool : Whether to save absorber information or not
        multi_plot_kwargs : dict : dictionary of arguments to pass to the multi_plot
                class

    Returns:
        none
    """

    #create initial slice
    mp = absorber_plotter(ray_filename=rays[0], **multi_plot_kwargs)
    mp.create_slice(cmap='cividis')

    #clear any annotations
    mp.slice.annotate_clear()

    for ray_fname in rays:
        #load in new ray
        mp.load_ray(ray_fname)

        #add ray and other annotations
        ray_num = get_ray_num(ray_fname)
        mp.add_annotations()
        mp.slice.annotate_title(f"ray {ray_num}")

        #create and save frame
        outfile = f"{out_dir}/mp{ray_num}.png"
        mp.create_multi_plot(outfile)

        if save_data:
            # save absorber data
            if mp.ice_table is not None:
                outfile=f"{out_dir}/ray{ray_num}_ice_absorbers{mp.num_ice}.h5"
                mp.ice_table.write(outfile, overwrite=True)

            if mp.spectacle_table is not None:
                outfile=f"{out_dir}/ray{ray_num}_spectacle_absorbers{mp.num_spectacle}.h5"
                mp.ice_table.write(outfile, overwrite=True)

        #close ray files and clear axes/annoations
        mp.ray.close()
        mp.fig.clear()
        mp.slice.annotate_clear()


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


if __name__ == '__main__':
    #take in arguments
    if len(argv) == 8:
        filename = argv[1]
        ray_dir = argv[2]
        ion_name = argv[3]
        out_dir= argv[4]
        use_bv = argv[5]
        frac = float(argv[6])
        cuts = argv[7]

    else:
        raise RuntimeError("Takes 7 arguments: Dataset_fname Ray_directory Ion_name Output_directory use_bv? frac cgm_cut? ")

    #check to see if should use bulk velocity
    if use_bv == 'True':
        use_bv=True
    else:
        use_bv=False

    #Define whether to use cuts. See Foggie's github. consistency files for
    #reason behind values
    #density in g/cm^3. Temp in K. radius in kpc
    if cuts == 'cgm':
        cut_list ="((obj[('gas', 'radius')].in_units('kpc') > 10) & \
                   (obj[('gas', 'radius')].in_units('kpc') < 200)) & \
                   ((obj[('gas', 'temperature')].in_units('K') > 1.5e4) | \
                   (obj[('gas', 'density')].in_units('g/cm**3') < 2e-26))"
    elif cuts == 'ism':
        cut_list =[["obj[('gas', 'radius')].in_units('kpc') < 10"],
                   ["obj[('gas', 'temperature')].in_units('K') < 1.5e4"],
                   ["obj[('gas', 'density')].in_units('g/cm**3') > 2e-26"]]
    else:
        cut_list=None
    main(filename, ray_dir, ion_name, out_dir, use_bv, frac, cut_list = cut_list)
