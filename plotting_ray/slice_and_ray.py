import yt
import trident as tri
import numpy as np

#filename=home+"data/DD0076/DD0076"
filename="/home/bb/Repo/yt/doc/source/quickstart/IsolatedGalaxy/galaxy0030/galaxy0030"
#load in dataset
ds = yt.load(filename)
center_ds = np.array([0.5, 0.5, 0.5])

ray_start= np.array([0.55, 0.35, 0.52], dtype=float)
ray_end = np.array([0.55, 0.65, 0.58], dtype=float)
nrth_vec=np.array([0, 0, 1])
line_list = ['H', 'Si', 'Mg II', 'C II 1335']


#find ray and normalize it
ray = ray_end - ray_start
ray = ray/np.linalg.norm(ray)

ray_center = (ray_end + ray_start)/2
#dif = np.array(ray_center) - np.array(center_ds)
#find ray perpendicular with z=0
norm_vec = np.cross(ray, nrth_vec)
#project center to plane of slice
center = np.dot(ray_center - center_ds, norm_vec) *norm_vec + center_ds
print(ray_center)
print(norm_vec)

print(center)
cut = yt.SlicePlot(ds, norm_vec, 'density', north_vector = nrth_vec, center=center)

ray_start =ds.arr(ray_start, "code_length")
ray_end = ds.arr(ray_end, "code_length")
tri_ray = tri.make_simple_ray(ds,
                              start_position = ray_start,
                              end_position = ray_end,
                              lines = line_list,
			      ftype = 'gas')

cut.annotate_ray(tri_ray)
cut.save("attempted_cut.png")
