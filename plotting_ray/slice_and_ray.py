import yt
import trident as tri

home="/mnt/home/boydbre1"
filename=home+"DD0076/DD0076"

#load in dataset
ds = yt.load(filename)

ray_start= np.array([0, 0, 0])
ray_end = np.array([1, 1, 1])

line_list = [‘H’, ‘Si’, ‘Mg II’, ‘C II 1335’]

#calc norm vec and center
center = (ray_end + ray_start)/2

#find ray and normalize it
ray = ray_end - ray_start
ray = np.linalg.norm(ray)

#find ray perpendicular with z=0
norm_vec = [ray[0], -ray[1], 0]


cut = yt.SlicePlot(ds, norm_vec, 'density', [0, 0, 1])

tri_ray = tri.make_simple_ray(ds,
                              start_position = ray_start,
                              end_position = ray_end,
                              lines = line_list,
                              ftype='gas')

cut.annotate_ray(ray)
cut.save("attempted_cut.png")
