import yt
import trident
import numpy as np
from mpi4py import MPI
from scipy.spatial.transform import Rotation
from os import makedirs
import seaborn as sns

def create_proj_frames(ds_fname,
                       center,
                       normal_vec,
                       offset=30,
                       fields=['density'],
                       color_maps=['magma'],
                       ccwh_gas=True,
                       weight=None,
                       num_frames=10,
                       width = 100,
                       out_dir="./"
                       ):
    """
    Creates projection plots at different rotations. Then frames can be combined
    to make a movie. Divides rotations at which projections are made evenly
    between number of proccess given.

    Parameters:
        ds_fname : str/path : path to the dataset
        center :list/array floats: coordiantes of the galaxy's center (in code_length)
        normal_vec :list/array floats: vector pointing perpendicular to galaxy's disk
        offset :float : angle to "tilt" the galaxy (in degrees)
        fields : list of str: Fields to make projections of
        color_maps : list of str: colormaps to use for each field. Must be same
                    size as fields and indices should match up.
        ccwh_gas : bool : Also Create projections of the density of cold, cool,
                    warm and hot gas. Cold gas <10^4 K. Cool 10^4 - 10^5 K.
                    Warm gas 10^5 - 10^6 K. Hot gas > 10^6 K.
        weight : str : field to use to make a weighted projection
        num_frames : int: number of frames to create for one full rotation
        width : float : width of the projection in kpc
        out_dir : str/path: path to where files will be saved

    Returns:
        none
    """

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
    lim_dict = dict(density=[5e-2, 1e4], #in msun/pc^2
                    temperature=[1e29, 1e31],
                    metallicity=[1e22, 1e24],
                    H_p0_number_density=[1e12, 1e24],
                    Si_p1_number_density=[1e10, 1e17],
                    Si_p2_number_density=[1e11, 1e16],
                    Si_p3_number_density=[1e11, 1e16],
                    C_p1_number_density=[1e11, 1e16],
                    C_p3_number_density=[1e11, 1e16],
                    O_p5_number_density=[1e11, 1e15],
                    cold=[1e-5, 1e-1],
                    cool=[1e-6, 1e-2],
                    warm=[1e-7, 1e-3],
                    hot=[1e-8, 1e-4])

    titles = dict(density="Density",temperature="Temperature", metallicity="Metallicity", H_p0_number_density="H I",
                    Si_p1_number_density="Si II",
                    Si_p2_number_density="Si III",
                    Si_p3_number_density="Si IV",
                    C_p1_number_density="C II",
                    C_p3_number_density="C IV",
                    O_p5_number_density="O VI")

    #load ds and construct sphere around galaxy
    ds = yt.load(ds_fname)
    trident.add_ion_fields(ds, ['Si II', 'Si III','Si IV', 'C II', 'C IV', 'O VI'])
    sph = ds.sphere(center, (width, 'kpc'))
    pad = int(np.ceil( np.log10(num_frames)))
    for i in my_rot_nums:
        #construct rotation vector and use to rotate
        rot_vec = normal_vec * rotations[i]
        rot = Rotation.from_rotvec(rot_vec)
        proj_vec = rot.apply(f_proj_vec)
        for fld, cmap in zip(fields, color_maps):
            prj = yt.OffAxisProjectionPlot(ds, proj_vec, fld,
                                           center=center, width=(width, 'kpc'),
                                           north_vector=normal_vec,
                                           weight_field=weight,
                                           data_source=sph)
            # set color bar and color map to be consistent on all proj
            lim_lb, lim_ub = lim_dict[fld]
            prj.set_unit('density', 'Msun/pc**2')
            prj.set_zlim(fld, lim_lb, lim_ub)
            prj.set_cmap(field=fld, cmap=cmap)
            prj.set_background_color(field=fld)
            prj.annotate_title(titles[fld])
            prj.annotate_scale()
            prj.hide_axes(draw_frame=True)
            prj.save(f"{out_dir}/{fld}/proj{i:0{pad}d}.png")

        if ccwh_gas:
            #make projections of dif temp gas
            names=['cold', 'cool', 'warm', 'hot']
            temps = [ [0, 1e4], [1e4, 1e5], [1e5, 1e6], [1e6, 1e10]]
            labels = ["Cold Gas Density $(T < 10^4 K)$", "Cool Gas Density $(10^4 < T < 10^5 K)$",
                      "Warm Gas Density $(10^5 < T < 10^6 K)$", "Hot Gas Density $(T > 10^6 K)$"]

            for temp, name, label in zip(temps, names, labels):
                reg = ds.cut_region(sph, [f"obj['temperature'] > {temp[0]}",
                                          f"obj['temperature'] < {temp[1]}"])
                prj = yt.OffAxisProjectionPlot(ds, proj_vec, 'density',
                                               center=center, width=(width, 'kpc'),
                                               north_vector=normal_vec,
                                               weight_field=weight,
                                               data_source=reg)

                lim_lb, lim_ub = lim_dict[name]
                lim_lb = ds.quan(lim_lb, 'g/cm**2').in_units('Msun/pc**2')
                lim_ub = ds.quan(lim_ub, 'g/cm**2').in_units('Msun/pc**2')
                prj.set_unit('density', 'Msun/pc**2')

                prj.set_zlim('density', lim_lb, lim_ub)
                prj.set_cmap(field='density', cmap='magma')
                prj.set_background_color('density')
                prj.annotate_title(label)
                prj.annotate_scale()
                prj.hide_axes(draw_frame=True)
                prj.save(f"{out_dir}/{name}_gas/proj{i:0{pad}d}.png")

if __name__ == '__main__':
    import sys
    sys.path.insert(0, '/mnt/home/boydbre1/Repo/CGM/')
    from multi_plot_movie.center_finder import find_center

    dsname = sys.argv[1]
    frms = int(sys.argv[2])
    offset=float(sys.argv[3])
    width=float(sys.argv[4])
    out_dir = sys.argv[5]

    #fields = ["C_p3_number_density", "O_p5_number_density"]
    #cmaps = ['magma', 'magma']
    #fields = ["density", "H_p0_number_density", "temperature", "metallicity"]
    #cmaps = ["magma", "magma", "thermal", "haline"]
    fields = ['density', 'H_p0_number_density', 'Si_p1_number_density', 'Si_p2_number_density', 'Si_p3_number_density', 'C_p1_number_density', 'C_p3_number_density', 'O_p5_number_density']
    h1_cmap = sns.blend_palette(("white", "#ababab", "#565656", "black",
                                  "#4575b4", "#984ea3", "#d73027",
                                  "darkorange", "#ffe34d"), as_cmap=True)
    cmaps = ['magma', h1_cmap, 'ice', 'ice', 'ice', 'cividis', 'cividis', 'haline']
    #cmaps= ["magma"] *len(fields) 
    c, n, r, bv = find_center(dsname)
    makedirs(out_dir, exist_ok=True)
    for f in fields:
        makedirs(f"{out_dir}/{f}", exist_ok=True)

    #names=['cold', 'cool', 'warm', 'hot']
    #for f in names:
    #    makedirs(f"{out_dir}/{f}_gas", exist_ok=True)
    create_proj_frames(dsname, c, n, ccwh_gas=False,fields=fields, color_maps=cmaps, num_frames=frms, width=width,offset=offset, out_dir=out_dir)
