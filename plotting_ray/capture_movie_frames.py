import plotter
from sys import argv

#setup conditions
line_list = ['H I', 'O VI', 'C IV']

filename = argv[1]

in_dir = argv[2]
i_name = argv[3]

out_dir= argv[4]

print(filename, in_dir, i_name, out_dir)
#now actual test
movie =plotter.movie_multi_plot(filename,in_dir, ion_name=i_name ,absorber_fields=line_list, out_dir=out_dir)
movie.create_movie()
print("--------------finished----------------")

