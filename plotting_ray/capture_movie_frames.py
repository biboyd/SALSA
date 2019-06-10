import plotter
from mpi4py import MPI
import numpy as np
from sys import argv

#init mpi
comm = MPI.COMM_WORLD

#setup conditions
line_list = []#['H I', 'O VI', 'C IV']

#take in arguments
if len(argv) == 5:
    filename = argv[1]
    in_dir = argv[2]
    i_name = argv[3]
    out_dir= argv[4]
else:
    raise RuntimeError("Takes 4 arguments: Dataset_fname Ray_directory Ion_name Output_directory")

#now actual test
movie =plotter.movie_multi_plot(filename, in_dir, ion_name=i_name,
				absorber_fields=line_list,
				out_dir=out_dir,
				wavelength_width=100,
				resolution=0.1)

#split up ray id numbers betweeen proccesors
rayrange = np.arange( len(movie.ray_files) )
rayrange_split = np.array_split(rayrange, comm.size)
my_range = rayrange_split[ comm.rank ]

#calc the number_density limits if root processor
if comm.rank == 0:
    #get middle ray to represent scale
    num_rays = len(movie.ray_files)
    middle_ray_file = movie.ray_files[ int(num_rays/2) ]
    mid_ray= yt.load( f"{movie.ray_dir}/{middle_ray_file}" )

    #get median num density
    num_density = np.array(mid_ray.data[ f"{movie.ion_p_name()}_number_density" ], dtype=np.float64)
    med = np.median(num_density)

    #estimate min max values to number dense plot. and markers positioning
    self.num_dense_min = 0.01*med
    self.num_dense_max = 1000*med

else:
    num_density_range = np.empty(0., dtype=np.float64)

comm.Barrier()
#broadcast the number density limits
comm.Bcast( [num_density_range, MPI.DOUBLE] )

#create movie frames
movie.create_movie(num_dense=num_density_range, ray_range = my_range)
print("-------------- {} finished----------------".format(comm.rank))

comm.Barrier()
if comm.rank==0:
    print("-"*76)
    print("All process finished :)")
    print("-"*76)
