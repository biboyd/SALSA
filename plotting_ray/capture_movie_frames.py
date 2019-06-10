import plotter
from mpi4py import MPI
import numpy
from sys import argv

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

if comm.rank == 0:
    num_density_range = fjdkfj

else:
    num_density_range = np.empty(0., dtype=np.float64)

comm.Barrier()
comm.Bcast( [num_density_range, MPI.DOUBLE] )

movie.create_movie(num_dense=num_density_range, ray_range = my_range)

print("-------------- {} finished----------------".format(comm.rank))

comm.Barrier()
if comm.rank==0:
    print("-"*76)
    print("All nodes finished :)")
    print("-"*76)
