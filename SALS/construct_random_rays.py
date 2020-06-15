import yt
import trident
import numpy as np

from yt.data_objects.static_output import \
    Dataset

from mpi4py import MPI
from sys import argv
from os import makedirs
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt

from SALS.utils.construct_rays import construct_rays
from SALS.utils.filter_definitions import radius_function, ion_p

def random_sightlines(ds_file, center, num_sightlines, max_impact_param, min_impact_param=0, length=200):
    """
    randomly sample impact parameter to get random sightlines from a given galaxy center

    Parameters:
        : ds_file: str or YT dataset
        center : arr: coordinates of the center of the galaxy
        num_sightlines : int : number of sightlines to return
        max_impact_param : float : maximum impact param to sample from in kpc
        min_impact_param : float : minimum impact param to sample from in kpc
        length : float : length of the sightline in kpc

    Returns :
        start_points : array : 2d array of the startpoints for
                    each sightline
        end_points : array : 2d array of the endpoints for
                    each sightline
    """

    #set file names and ion name
    if isinstance(ds_file, str):
        ds = yt.load(ds_file)
    elif isinstance(ds_file, Dataset):
        ds = ds_file

    length = ds.quan(length, 'kpc').in_units('code_length')
    length = length.value
    min_impact_param = ds.quan(min_impact_param, 'kpc').in_units('code_length')
    max_impact_param = ds.quan(max_impact_param, 'kpc').in_units('code_length')

    #randomly select angle and distance from center of gal
    #take sqrt so that impact param is uniform in projected area space
    impact_param = np.sqrt(np.random.uniform(min_impact_param.value**2, max_impact_param.value**2, num_sightlines))

    #theta represents polar angle. phi represents azimuthal
    theta = np.random.uniform(0, np.pi, num_sightlines)
    phi = np.random.uniform(0, 2*np.pi, num_sightlines)

    #construct vector from gal_center to sightline midpoint
    rad_vec= np.empty((num_sightlines, 3))
    rad_vec[:, 0] = impact_param*np.cos(phi)*np.sin(theta)
    rad_vec[:, 1] = impact_param*np.sin(phi)*np.sin(theta)
    rad_vec[:, 2] = impact_param*np.cos(theta)

    #define vector along sightline (perpendicular to radial vector)
    perp_vec = np.empty_like(rad_vec)
    perp_vec[:, 0] = rad_vec[:, 1]
    perp_vec[:, 1] = -1*rad_vec[:, 0]
    perp_vec[:, 2] = 0.

    #randomly rotate perp_vec around rad vec
    alpha=np.random.uniform(0., 2*np.pi, num_sightlines)
    for i in range(num_sightlines):
        #normalize perpendicular vector
        perp_vec[i, :] =  perp_vec[i, :]/np.sqrt(perp_vec[i, 0]**2 + perp_vec[i, 1]**2)
        #rotate around radial vector by random amount
        rot_vec = alpha[i] *rad_vec[i, :]/np.linalg.norm(rad_vec[i, :])
        rot = Rotation.from_rotvec(rot_vec)
        perp_vec[i, :] = rot.apply(perp_vec[i, :])


    #shift to be centered at galaxy
    sightline_centers = rad_vec +np.array(center)

    #find ending and start points for each sightline
    end_point = sightline_centers + length/2 *perp_vec
    start_point = sightline_centers - length/2 *perp_vec

    return start_point, end_point, impact_param

def generate_lrays(ds, center,
                n_rays, max_impact_param,
                min_impact_param=0.,
                length=200,
                fld_params={},
                ion_list=['H I', 'C IV', 'O VI'],
                fields=['density','temperature'],
                out_dir='./'):
    """

    Parameters
    ----------
        :ds_file : str or YT Dataset
            path to dataset or already loaded, YT dataset

        :center : arr
            coordinates of the center of the galaxy

        :n_rays : int
            number of light rays to construct

        :max_impact_param : float
            maximum impact param to sample from in kpc

        :min_impact_param : float
            minimum impact param to sample from in kpc

        :length : float
            length of the sightline in kpc

        :ion_list : list
            ions to add to lightray

        :fields : list
            fields to add to lightray

        :out_dir : string
            path to where ray files will be written

    Returns
    -------
        none
    """

    comm = MPI.COMM_WORLD
    #collect sightlines
    if comm.rank == 0:
        start_pnts, end_pnts, imp_param = random_sightlines(ds, center,
                                                 n_rays,
                                                 max_impact_param,
                                                 min_impact_param=min_impact_param,
                                                 length=length)
        imp_param = ds.arr(imp_param, 'code_length').in_units('kpc')
        np.save(f"{out_dir}/impact_parameter.npy", imp_param)

    else:
        start_pnts= np.empty(n_rays, dtype=np.float64)
        end_pnts= np.empty(n_rays, dtype=np.float64)
        #imp_param= np.empty(n_rays, dtype=np.float64)

    # share sightline points
    comm.Barrier()
    comm.Bcast([start_pnts, MPI.DOUBLE])
    comm.Bcast([end_pnts, MPI.DOUBLE])
    #comm.Bcast([imp_param, MPI.DOUBLE])

    #add center to field parameters
    fld_params['center']=center

    #add ion fields to dataset if not already there
    trident.add_ion_fields(ds, ions=ion_list)

    construct_rays(ds, start_pnts, end_pnts,
                   fld_params=fld_params,
                   line_list=ion_list,
                   other_fields=fields,
                   out_dir=out_dir)
