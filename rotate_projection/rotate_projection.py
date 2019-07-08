import yt
import numpy as np
from mpi4py import MPI
from scipy.spatial.transform import Rotation

def create_proj_frames(ds_fname,
                       center,
                       normal_vec,
                       offset=30,
                       field='density',
                       weight=None,
                       num_frames=10,
                       out_dir="./"
                       ):

    comm = MPI.COMM_WORLD

    if comm.rank == 0:
        #get in plane of galaxy vector to rot offset
        inplane_vec = np.array( [norm_vec[2], 0, -norm_vec[0]] )

        #check not zero vector
        if inplane_vec[0] == 0 and inplane_vec[2] == 0:
            #switch to a non zero vector
            inplane_vec = np.array( [0, norm_vec[2], -norm_vec[1]] )

        inplane_vec = inplane_vec/np.linalg.norm(inplane_vec)
        #create rotation vector then rotate norm vec to get our axis of rot
        rot_vec = inplane_vec * np.deg2rad(offset)
        rot = Rotation.from_rotvec(rot_vector)
        axis_rot_vec = rot.apply(normal_vec)
        axis_rot_vec = np.float64(axis_rot_vec)

    else:
        axis_rot_vec = np.zeros(3, dtype=np.float64)
        inplane_vec = np.zeros(3, dtype=np.float64)

    comm.Barrier()
    comm.Bcast([axis_rot_vec, MPI.DOUBLE])
    comm.Bcast([inplane_vec, MPI.DOUBLE])
    rotatations = np.linspace(0, 2*np.pi, num_frames)
    my_rot_nums = np.arange(num_frames)
    #split ray numbers then take a portion based on rank
    split_rot_nums = np.array_split(my_rot_nums, comm.size)
    my_rot_nums = split_rot_nums[comm.rank]

    #load ds and construct sphere around galaxy
    ds = yt.load(ds_fname)
    sph = ds.sphere(center, (50, 'kpc'))

    pad = np.ceil( np.log10(num_frames) )
    for i in my_rot_nums:
        #construct rotation vector and use to rotate
        rot_vec = axis_rot_vec * rotations[i]
        rot = Rotation.from_rotvec(rot_vec)
        proj_vec = rot.apply(inplane_vec)

        prj = yt.OffAxisProjectionPlot(ds, proj_vec, field,
                                       center=center, width=(75, 'kpc'),
                                       north_vector=normal_vec,
                                       weight_field=weight,
                                       data_source=sph)
        prj.save(f"{out_dir}/proj{i:0{pad}d}.png")

if __name__ == '__main__':
    dsname = argv[1]

    c, n, r, bv = find_center(dsname)
    create_proj_frames('x', 1, 2)
