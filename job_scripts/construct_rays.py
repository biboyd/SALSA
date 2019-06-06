from sys import argv, path
path.append("/mnt/home/boydbre1/Repo/CGM/plotting_ray/")
import plotter

#setup conditions
line_list = ['H I', 'O VI', 'C IV']

filename = argv[1]

n_rays=int(argv[2])
ray_length=int(argv[3])
out_dir = argv[4] 

#now actual test
plotter.construct_rays(filename, line_list, 
			n_rays=n_rays, 
			length=ray_length, 
			dist_from_center=ray_length/2,
			out_dir=out_dir)

print("finished")
