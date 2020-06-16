import numpy as np
import yt
import trident
import pandas as pd

from salsa.absorber_extractor import absorber_extractor
from salsa.utils.collect_files import collect_files, check_rays
from salsa.utils.filter_definitions import ion_p_num
from salsa.generate_light_rays import generate_lrays
from mpi4py import MPI

from yt.data_objects.static_output import \
    Dataset

def generate_catalog(ds_file, n_rays,
                     ray_directory,
                     ion_list,
                     method="ice",
                     center=None,
                     impact_param_lims=(0, 200),
                     ray_length=200,
                     field_parameters={},
                     fields=[],
                     cut_region_filters=[],
                     extractor_kwargs={},
                     units_dict={}):

    """
    Generates a catalog of absorber properties from a given number of lightrays
    through a dataset. Will look if lrays have already been made, otherwise will
    create them by uniform randomly sampling impact parameter. Uses OpenMPI to
    split light ray creation and absorber extraction among processors.

    Parameters
    ----------
    ds_file: str or dataset
        either path to a dataset or a loaded dataset

    n_rays: int
        number of rays to sample

    ray_directory: str
        path to directory where rays loaded from or will be saved if they haven't
        been constructed

    ion_list: list str
        list of ions to find absorbers from.

    method: "ice" or "spectacle", optional
        Choose which method to use to extract absorbers.

    center: list or array, optional
        The center of the galaxy in units 'code_length'. If None, defaults to
        domain_center.
        Default: None

    impact_param_lims: tuple or list, optional
        The range on which to sample impact parameter when constructing lightrays
        Default: (0, 200)

    ray_length: float, optional
        The length of each light ray in units kpc.
        Default: 200

    field_parameters: dict, optional
        The parameters that will be passed to trident during ray construction.
        This can be something like "bulk_velocity" so that a radial velocity can
        be saved.

    fields: list of str
        YT fields to add to lightrays. Will be included in catalog if "ice" method
        is selected

    cut_region_filters: list of strings, optional
        a list of filters defined by the way you use Cut Regions in YT
        Default: None

    extractor_kwargs: dict, optional
        Additional key word arguments to pass to the absorber_extractor to
        modify default extraction parameters
        Default: {}

    units_dict: dict, optional
        dictionary of units to use for the fields when extracting properties
        (only relevant for 'ice' method)
        Default: None
    Returns
    -------
    full_catalog: pandas.DataFrame
        pandas dataframe containing all of the absorbers extracted from all
        the lightrays
    """
    comm = MPI.COMM_WORLD

    # check ds_file is need to load
    if isinstance(ds_file, str):
        ds = yt.load(ds_file)
    elif isinstance(ds_file, Dataset):
        ds = ds_file


    # add ion number density to fields to check
    check_fields = fields.copy()
    for i in ion_list:
        check_fields.append(ion_p_num(i))

    #check if rays already made
    if comm.rank == 0:
        check =check_rays(ray_directory, n_rays, check_fields)
        ray_bool= np.array([check], dtype=int)
    else:
        ray_bool = np.array([False], dtype=int)

    # share if rays made already or not
    comm.Barrier()
    comm.Bcast([ray_bool, MPI.INT])

    #generate rays randomly
    if not ray_bool[0]:
        #set a center
        if center is None:
            center=ds.domain_center

        #construct random rays in ray_directory
        generate_lrays(ds, center,
                    n_rays, impact_param_lims[1],
                    min_impact_param=impact_param_lims[0],
                    length=ray_length,
                    fld_params=field_parameters,
                    ion_list=ion_list,
                    fields=fields,
                    out_dir=ray_directory)

    #Extract Absorbers

    #collect and split up ray files
    ray_files = np.array(collect_files(ray_directory, key_words=['ray']), dtype=str)
    ray_files_split = np.array_split(ray_files, comm.size)
    my_rays = ray_files_split[ comm.rank ]

    #add directory path to rays
    my_ray_files=[ ray_directory+'/'+r for r in my_rays ]

    #create catalog for each ion
    df_list=[]
    for ion in ion_list:
        # setup absorber extractor
        abs_ext = absorber_extractor(ds, my_ray_files[0], ion_name=ion,
                                     cut_region_filters=cut_region_filters,
                                     **extractor_kwargs)

        # get catalogs
        my_df = get_catalog(abs_ext, my_ray_files, method, fields=fields, units_dict=units_dict)
        if my_df is not None:
            df_list.append(my_df)

    my_catalog= pd.concat(df_list)
    comm.Barrier()

    #gather all catalogs and creae one large
    all_dfs= comm.allgather(my_catalog)
    full_catalog = pd.concat(all_dfs)

    return full_catalog



def get_catalog(abs_extractor, ray_list, method, fields=None, units_dict=None):
    """
    Create catalog of the given rays usin absorber extractor

    Parameters
    ----------
    abs_extractor: SALS.absorber_extractor
        Absorber Extractor object that will be used to extract absoprtion feat.
        for the catalog

    ray_list: list of str or trident.ray objects
        List of ray objects or list of trident rays whose absorbers will be
        extracted
    method: str
        Either 'ice' or 'spectacle'. specifies which method is used to extract
        absorbers

    fields: list str, optional
        Fields to extract/add to catalog if using 'ice' method.
        Defaults=None

    units_dict: dict
        dictionary containing what units to use for each field

    Returns
    -------
    full_df: pandas.DataFrame
        Catalog of absorber properties in a pandas dataframe.
    """
    df_list=[]

    if method == 'ice':
        for ray in ray_list:
            #load new ray and extract absorbers
            abs_extractor.load_ray(ray)
            df = abs_extractor.get_ice_absorbers(fields=fields, user_unit_dict=units_dict)

            if df is not None:
                # add ray index
                ray_num = get_ray_num(ray)
                start = 65 # Ascii number for 'A'
                for i in range(abs_extractor.num_ice):
                    df.loc[i,'absorber_index'] = f"{ray_num}{chr(start+i)}"
                df_list.append(df)

    elif method == 'spectacle':
        for ray in ray_list:
            abs_extractor.load_ray(ray)
            df = abs_extractor.get_spectacle_absorbers()

            #add ray index
            if df is not None:
                ray_num = get_ray_num(ray)
                start = 65 # Ascii number for 'A'
                for i in range(abs_extractor.num_spectacle):
                    df.loc[i,'absorber_index'] = f"{ray_num}{chr(start+i)}"
                df_list.append(df)

    else:
        raise RuntimeError(f"method={method} is not valid. method must be 'ice' or 'spectacle'.")

    if len(df_list) > 0:
        full_df = pd.concat(df_list)
    else:
        full_df = None
    return full_df


def get_ray_num(file_path):
    """
    extract the ray's number from it's file name by removing 'ray' and '.h5' as
    well as preceding path
    """
    filename = file_path.split('/')[-1]
    num = filename[3:-3]
    return num
