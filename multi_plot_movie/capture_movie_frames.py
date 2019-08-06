from mpi4py import MPI
import numpy as np
from sys import argv
import yt
import trident
from multi_plot import multi_plot
from center_finder import find_center
from os import makedirs, listdir

def main():
    #init mpi
    comm = MPI.COMM_WORLD

    #setup conditions
    line_list = []#['H I', 'O VI', 'C IV']

    #take in arguments
    if len(argv) == 7:
        filename = argv[1]
        ray_dir = argv[2]
        i_name = argv[3]
        out_dir= argv[4]
        use_bv = argv[5]
        sigma = float(argv[6])
    else:
        raise RuntimeError("Takes 5 arguments: Dataset_fname Ray_directory Ion_name Output_directory use_bv?")

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
    comm.Bcast([bulk_vel, MPI.DOUBLE])
    comm.Barrier()

    #check to see if should use bulk velocity
    if use_bv != 'True':
        bulk_vel=None

    #set up multiplot settings
    mp_kwargs = dict(ds_filename=filename, ion_name=i_name,
    				    absorber_fields=line_list,
                                    center_gal = center,
                                    north_vector = nvec,
                                    redshift=rshift[0],
                                    bulk_velocity=bulk_vel,
                                    use_spectacle= True,
                                    plot_spectacle=True,
                                    sigma_smooth= sigma,
                                    wavelength_width=20,
                                    resolution=0.1)


    #collect only ray files in ray_dir
    ray_files=[]
    for f in listdir(ray_dir):
        #check hdf5 file
        if (f[-3:] == ".h5"):
            full_name ="/".join((ray_dir, f))
            ray_files.append(full_name)

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
                      'num_dense_max': num_density_range[1],
                      'markers_nd_pos': num_density_range[0]*5})

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
                  multi_plot_kwargs={}):
    """
    creates a movie by combining all the plots made from the ray in ray_dir

    Parameters:
        rays
        num_dense
        slice_width
        slice_height
    """

    #create initial slice
    mp = multi_plot(ray_filename=rays[0], **multi_plot_kwargs)
    mp.create_slice()

    #clear annotations
    mp.slice.annotate_clear()

    for ray_fname in rays:
        #load in new ray
        mp.ray = yt.load(ray_fname)

        #add ray and other annotations
        ray_num = get_ray_num(ray_fname)
        mp.add_annotations()
        mp.slice.annotate_title(f"ray {ray_num} sigma={mp.sigma_smooth}")

        #create and save frame
        outfile = f"{out_dir}/mp{ray_num}.png"
        mp.create_multi_plot(outfile)

        #close ray files and clear axes/annoations
        mp.ray.close()

        mp.fig.clear()
        mp.slice.annotate_clear()


def get_ray_num(file_path):
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
    main()
