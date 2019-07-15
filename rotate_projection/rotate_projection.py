import sys
sys.path.insert(0, '/mnt/home/boydbre1/Repo/CGM/cosmo_analysis/')
import yt
import trident
import numpy as np
from mpi4py import MPI
from scipy.spatial.transform import Rotation
from os import makedirs
from center_finder import find_center

def create_proj_frames(ds_fname,
                       center,
                       normal_vec,
                       offset=30,
                       fields=['density'],
                       color_maps=['magma'],
                       weight=None,
                       num_frames=10,
                       out_dir="./"
                       ):

    comm = MPI.COMM_WORLD
    normal_vec = np.array(normal_vec)
    if comm.rank == 0:
        #get in plane of galaxy vector to rot offset
        inplane_vec = np.array( [normal_vec[2], 0, -normal_vec[0]] )

        #check not zero vector
        if inplane_vec[0] == 0 and inplane_vec[2] == 0:
            #switch to a non zero vector
            inplane_vec = np.array( [0, normal_vec[2], -normal_vec[1]] )

        inplane_vec = inplane_vec/np.linalg.norm(inplane_vec)
        #create rotation vector then rotate norm vec to get our axis of rot
        rot_vec = inplane_vec * np.deg2rad(offset)
        rot = Rotation.from_rotvec(rot_vec)
        off_axis_vec = rot.apply(normal_vec)
        off_axis_vec = np.float64(off_axis_vec)

        f_proj_vec = np.cross(off_axis_vec, inplane_vec)
        f_proj_vec = f_proj_vec/np.linalg.norm(f_proj_vec)

    else:
        f_proj_vec = np.zeros(3, dtype=np.float64)

    comm.Barrier()
    comm.Bcast([f_proj_vec, MPI.DOUBLE])
    interval = 2*np.pi/num_frames
    rotations = np.arange(0, 2*np.pi, interval)
    my_rot_nums = np.arange(num_frames)
    #split ray numbers then take a portion based on rank
    split_rot_nums = np.array_split(my_rot_nums, comm.size)
    my_rot_nums = split_rot_nums[comm.rank]

    #limits dictionary to fix colobar scale
    lim_dict = dict(density=[1e-6, 1e0],
                    temperature=[1e29, 1e31],
                    metallicity=[1e22, 1e24],
                    H_p0_number_density=[1e15, 1e22],
                    C_p3_number_density=[1e10, 1e16],
                    O_p5_number_density=[1e14, 1e16],
                    cold=[1e-5, 1e-1],
                    cool=[1e-6, 1e-2],
                    warm=[1e-7, 1e-3],
                    hot=[1e-8, 1e-4])
    #load ds and construct sphere around galaxy
    ds = yt.load(ds_fname)
    trident.add_ion_fields(ds, ['C IV', 'O VI'])
    sph = ds.sphere(center, (100, 'kpc'))
    pad = int(np.ceil( np.log10(num_frames)))
    for i in my_rot_nums:
        #construct rotation vector and use to rotate
        rot_vec = normal_vec * rotations[i]
        rot = Rotation.from_rotvec(rot_vec)
        proj_vec = rot.apply(f_proj_vec)
        for fld, cmap in zip(fields, color_maps):
            break
            prj = yt.OffAxisProjectionPlot(ds, proj_vec, fld,
                                           center=center, width=(100, 'kpc'),
                                           north_vector=normal_vec,
                                           weight_field=weight,
                                           data_source=sph)
            # set color bar and color map to be consistent on all proj
            lim_lb, lim_ub = lim_dict[fld]
            prj.set_zlim(fld, lim_lb, lim_ub)
            prj.set_cmap(field=fld, cmap=cmap)
            prj.save(f"{out_dir}/{fld}/proj{i:0{pad}d}.png")
        
        #make projections of dif temp gas
        names=['cold', 'cool', 'warm', 'hot']
        temps = [ [0, 1e4], [1e4, 1e5], [1e5, 1e6], [1e6, 1e10]]
        labels = ["Cold Gas Density $(T < 10^4 K)$", "Cool Gas Density $(10^4 < T < 10^5 K)$",
                  "Warm Gas Density $(10^5 < T < 10^6 K)$", "Hot Gas Density $(T > 10^6 K)$"]

        for temp, name, label in zip(temps, names, labels):
            reg = ds.cut_region(sph, [f"obj['temperature'] > {temp[0]}", 
                                      f"obj['temperature'] < {temp[1]}"])
            prj = yt.OffAxisProjectionPlot(ds, proj_vec, 'density',
                                           center=center, width=(100, 'kpc'),
                                           north_vector=normal_vec,
                                           weight_field=weight,
                                           data_source=reg)
            
            lim_lb, lim_ub = lim_dict[name]
            prj.set_zlim('density', lim_lb, lim_ub)
            prj.set_cmap(field='density', cmap='magma')
            prj.set_background_color('density')
            prj.annotate_title(label)
            prj.annotate_scale()
            prj.hide_axes(draw_frame=True)
            prj.save(f"{out_dir}/{name}_gas/proj{i:0{pad}d}.png")

if __name__ == '__main__':
    dsname = sys.argv[1]
    frms = int(sys.argv[2])
    out_dir = sys.argv[3]

    #fields = ["C_p3_number_density", "O_p5_number_density"]
    #cmaps = ['magma', 'magma']
    fields = ["density", "H_p0_number_density", "temperature", "metallicity"]
    cmaps = ["magma", "magma", "thermal", "haline"]
    c, n, r, bv = find_center(dsname)
    makedirs(out_dir, exist_ok=True)
    for f in fields:
        makedirs(f"{out_dir}/{f}", exist_ok=True)
        
    names=['cold', 'cool', 'warm', 'hot']
    for f in names:
        makedirs(f"{out_dir}/{f}_gas", exist_ok=True)
    create_proj_frames(dsname, c, n, fields=fields, color_maps=cmaps, num_frames=frms, out_dir=out_dir)
