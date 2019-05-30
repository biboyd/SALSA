import plotter
from sys import argv

#setup conditions
line_list = ['H I', 'O VI', 'C IV']

filename = argv[1]

n_rays=int(argv[2])
out_dir = argv[3] 

#now actual test
plotter.construct_rays(filename, line_list, n_rays=n_rays, length=100, out_dir=out_dir)
print("finished")
