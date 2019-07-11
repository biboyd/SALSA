#creates slices of the galaxy as well as a density proj
import yt
import numpy as np
from center_finder import find_center
from sys import argv
from os import makedirs
from scipy.spatial.transform import Rotation

def main(dataset,
         length=200,
         norm_vector=[0, 0, 1],
         angle=0,
         bulk_vel = 0,
         center=None,
         out_dir='./info_slices'):

    #start yt parallelization
    yt.enable_parallelism()
    ds = yt.load(dataset)

    #create directory if doesn't exist
    makedirs(out_dir, exist_ok=True)



    #convert lengths to code_length and angle to radians
    length = ds.quan(length, 'kpc').in_units('code_length')
    norm_vector = np.array(norm_vector)/np.linalg.norm(norm_vector)
    #create the unit vector for rays
    ray_unit = np.array( [norm_vector[2], 0, -norm_vector[0]] )

    #check not zero vector
    if ray_unit[0] == 0 and ray_unit[2] == 0:
        #switch to a non zero vector
        ray_unit = np.array( [0, norm_vector[2], -norm_vector[1]] )

    ray_unit = ray_unit/np.linalg.norm(ray_unit)

    #rotate ray unit
    if angle != 0:
        angle = np.deg2rad(angle)
        rot_vector = norm_vector*angle
        rot = Rotation.from_rotvec(rot_vector)
        ray_unit = rot.apply(ray_unit)

    if center is None:
        center = ds.domain_center.in_units('code_length')
    else:
        center = ds.arr(center, 'code_length')


    #plot slices for density, temp and metallicity to compare with multi plot
    fields = ['density', 'temperature', 'metallicity', 'velocity_magnitude', 'vel_bv']
    color_maps = ['magma', 'thermal', 'haline', 'viridis', 'magma']
     
    #construct sphere to make slices/projections from
    sph = ds.sphere(center, length)
    for fld, cmap in yt.parallel_objects(zip(fields, color_maps)):
        #construct vector normal to slicing plane
        slc_norm = np.cross(ray_unit, norm_vector)
        if fld == 'vel_bv':
            #create vector field, take in account the galaxy's bulk velocity
            sph.set_field_parameter('bulk_velocity', ds.arr(bulk_vel, 'km/s'))
            slc_bv = yt.OffAxisSlicePlot(ds, slc_norm, 'density',
                                   north_vector = norm_vector,
                                   center=center, width=length,
                                   data_source=sph)
            slc_bv.set_axes_unit('kpc')
            slc_bv.set_cmap(field='density', cmap=cmap)
            slc_bv.set_background_color('density')

            #over plot new velocities
            slc_bv.annotate_quiver('cutting_plane_velocity_x', 'cutting_plane_velocity_y',
                                factor=24, plot_args={'color':'white'})
            slc_bv.annotate_title("Velocity Field in galaxy's reference frame")
            slc_bv.save(f"{out_dir}/velocity_field_bv.png")

        else:
            #crete slice of corresponding field
            slc = yt.OffAxisSlicePlot(ds, slc_norm, fld,
                                   north_vector = norm_vector,
                                   center=center, width=length,
                                   data_source=sph)
            slc.set_axes_unit('kpc')
            slc.set_cmap(field=fld, cmap=cmap)
            slc.set_background_color(fld)
            if fld == 'velocity_magnitude':
                slc.set_unit('velocity_magnitude', 'km/s')
            slc.save(f"{out_dir}/{fld}_slice.png")

            if fld == 'density':
                #overplot velocities
                slc.annotate_quiver('cutting_plane_velocity_x', 'cutting_plane_velocity_y',
                                    factor=24, plot_args={'color':'white'})
                slc.annotate_title("Velocity Field in observors reference frame")
                slc.save(f"{out_dir}/velocity_field_no_bv.png")

    #create projection plot
    prj = yt.OffAxisProjectionPlot(ds, norm_vector, 'density',
                                   center=center, width=length,
                                   data_source=sph)
    prj.set_axes_unit('kpc')
    prj.save(f"{out_dir}/density_projection.png")

if __name__ == '__main__':
    ds_fname = argv[1]
    length=int(argv[2])
    out_dir = argv[3]
    c, nv, r, bv = find_center(ds_fname)
    main(ds_fname, length=length,norm_vector=nv, bulk_vel=bv, center=c, out_dir=out_dir)
