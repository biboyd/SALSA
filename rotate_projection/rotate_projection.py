import sys
sys.path.insert(0, '/mnt/home/boydbre1/Repo/CGM/cosmo_analysis/')
import yt
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
    if comm.rank >= 0:
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

    #comm.Barrier()
    #comm.Bcast([f_proj_vec, MPI.DOUBLE])
    rotations = np.linspace(0, 2*np.pi, num_frames)
    my_rot_nums = np.arange(num_frames)
    #split ray numbers then take a portion based on rank
    split_rot_nums = np.array_split(my_rot_nums, comm.size)
    my_rot_nums = split_rot_nums[comm.rank]

    #load ds and construct sphere around galaxy
    ds = yt.load(ds_fname)
    sph = ds.sphere(center, (100, 'kpc'))
    pad = int(np.ceil( np.log10(num_frames)))
    for i in my_rot_nums:
        #construct rotation vector and use to rotate
        rot_vec = normal_vec * rotations[i]
        rot = Rotation.from_rotvec(rot_vec)
        proj_vec = rot.apply(f_proj_vec)
        for fld, cmap in zip(fields, color_maps):
            prj = yt.OffAxisProjectionPlot(ds, proj_vec, fld,
                                           center=center, width=(100, 'kpc'),
                                           north_vector=normal_vec,
                                           weight_field=weight,
                                           data_source=sph)
            prj.set_cmap(field=fld, cmap=cmap)
            prj.save(f"{out_dir}/{fld}/proj{i:0{pad}d}.png")

if __name__ == '__main__':
    dsname = sys.argv[1]
    frms = int(sys.argv[2])
    out_dir = sys.argv[3]

    fields = ["density", "temperature", "metallicity", "H_p0_number_density"]
    cmaps = ["magma", "thermal", "haline", "magma"]
    c, n, r, bv = find_center(dsname)
    makedirs(out_dir, exist_ok=True)
    for f in fields:
        makedirs(f"{out_dir}/{f}", exist_ok=True)
    create_proj_frames(dsname, c, n, fields=fields, color_maps=cmaps, num_frames=frms, out_dir=out_dir)
