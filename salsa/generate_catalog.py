import numpy as np
import yt
from astropy.table import vstack

from salsa.absorber_extractor import SPICEAbsorberExtractor
from salsa.utils.collect_files import collect_files, check_rays
from salsa.utils.functions import ion_p_num
from salsa.generate_light_rays import generate_lrays
from mpi4py import MPI

from yt.data_objects.static_output import \
    Dataset

def generate_catalog(ds_file, n_rays,
                     ray_directory,
                     ion_list,
                     method="spice",
                     center=None,
                     impact_param_lims=(0, 200),
                     ray_length=200,
                     field_parameters={},
                     fields=[],
                     ftype='gas',
                     extractor_kwargs={},
                     units_dict={}):

    """
    Generates a catalog of absorber properties from a given number of lightrays
    through a dataset. Will look if lrays have already been made, otherwise will
    create them by uniform randomly sampling impact parameter. Uses OpenMPI to
    split up light ray creation and absorber extraction among processors.

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

    method: "spice", optional
        Choose which method to use to extract absorbers.
        Currently only "spice" is supported.

    center: list or array, optional
        The center of the galaxy in units 'code_length'. If None, defaults to
        domain_center.
        Default: None

    impact_param_lims: tuple or list, optional
        The range on which to sample impact parameter when constructing lightrays.
        If no units are associated with the numbers (either astropy or unyt),
        the code will assume kpc.
        Default: (0, 200)

    ray_length: float, optional
        The length of each light ray. If no units are associated (either astropy
        or unyt) the code will assume kpc.
        Default: 200

    field_parameters: dict, optional
        The parameters that will be passed to trident during ray construction.
        This can be something like "bulk_velocity" so that a radial velocity can
        be saved.

    fields: list of str
        YT fields to add to lightrays. Will be included in catalog if "spice" method
        is selected

    ftype : str
        The field to be passed to trident that ion fields will be added to, i.e.
        ``('gas', 'H_p0_number_density')``. ``'gas'`` should work for most grid-based
        simulations. For particle-based simulations this will not work and needs
        to be changed. ``'PartType0'`` often works though it varies.
        See ``trident.add_ion_fields()`` for more information

    extractor_kwargs: dict, optional
        Additional key word arguments to pass to the absorber_extractor to
        modify default extraction parameters. Either a single dict that will be
        passed for each ion. Or a dict of ions pointing toward individual extractor
        kwargs. Examples:

        ``extractor_kwargs={'H I':{'absorber_min':14}, 'C IV':{'absorber_min':13}, 'O VI':{}}``
        ``extractor_kwargs={'absorber_min':13.5}``

        The first will set different absober mins for each ion, with O VI taking
        default as specified by ``salsa.utils.defaults.default_cloud_dict``. The
        second example will set the minimum absorber as 13.5 for every ion.
        **NOTE** you cannot mix the two formats. If one ion is specified then
        all ions must be specified (see 'O VI' included even though it's dictionary is empty)

        Default: {}

    units_dict: dict, optional
        dictionary of astropy units to use for the fields when extracting properties
        (only relevant for 'spice' method)
        Default: None

    Returns
    -------
    full_catalog: astropy QTable
        Table containing all of the absorbers extracted from all
        the lightrays. If no absorbers are found, None is returned
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

    # check if rays already made
    check = check_rays(ray_directory, n_rays, check_fields)
    my_ray_bool = np.array([check], dtype=int)
    ray_bool = np.array([0], dtype=int)

    # share if rays made already or not
    comm.Barrier()
    comm.Allreduce([my_ray_bool, MPI.INT],[ray_bool, MPI.INT], op=MPI.LAND)

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
                    ftype=ftype,
                    ray_directory=ray_directory)

    comm.Barrier()
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
        #check if extractor kwargs has ion specific information
        if ion in extractor_kwargs.keys():
            curr_kwargs = extractor_kwargs[ion]
        else:
            curr_kwargs = extractor_kwargs.copy()

        # setup absorber extractor
        if method == "spice":
            abs_ext = SPICEAbsorberExtractor(ds, 
                                            ion_name=ion,
                                            **curr_kwargs)
        else:
            raise RuntimeError(f"Method {method} is not recognized.")

        # get catalogs
        my_df = abs_ext.get_all_absorbers(my_ray_files,
                                          fields=fields,
                                          units_dict=units_dict)
        if my_df is not None:
            df_list.append(my_df)

    # Return Nonetype if no absorbers found
    if df_list == []:
        my_catalog = None
    else:
        my_catalog= vstack(df_list)
    comm.Barrier()

    #gather all catalogs and creae one large
    all_dfs= comm.allgather(my_catalog)

    # check if any absorbers were found
    if all(v is None for v in all_dfs):
        full_catalog = None
    else:
        full_catalog = vstack(all_dfs)

    return full_catalog
