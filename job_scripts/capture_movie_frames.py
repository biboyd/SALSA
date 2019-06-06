from sys import argv, path
path.append("/mnt/home/boydbre1/Repo/CGM/plotting_ray/")
import plotter

#setup conditions
line_list = ['H I', 'O VI', 'C IV']

#take in arguments
if len(argv) == 5:
    filename = argv[1]
    in_dir = argv[2]
    i_name = argv[3]
    out_dir= argv[4]
else:
    raise RuntimeError("Takes 4 arguments: Dataset_fname Ray_directory Ion_name Output_directory")

#now actual test
movie =plotter.movie_multi_plot(filename,in_dir, ion_name=i_name ,absorber_fields=line_list, out_dir=out_dir)
movie.create_movie()
print("--------------finished----------------")
