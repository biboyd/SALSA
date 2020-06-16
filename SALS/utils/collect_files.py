from sys import argv
from os import listdir
from astropy.table import QTable, vstack

import pandas as pd

def collect_files(directory, file_ext='.h5', key_words=[], black_list=[]):
    """
    Finds and returns files in a given directory with given extension. Optional
    key_word and black_list requirements

    Parameters
    ----------

    : directory : str
        Directory in which to search for files

    : file_ext : str
        The file extension to look for (ie '.py', '.h5')

    : key_words : list of str
        Key words that must in the filename to be collected

    : black_list : list of str
        list of files to exclude from collection

    Returns
    -------

    : files : list of str
        List of filenames in the directory that pass all requirements for file
        extension, keywords, black_list.

    """

    all_files = listdir(directory)

    # get np files only
    files=[]
    for f in all_files:
        #check has file extension and not black listed
        if file_ext in f and check_file(f, key_words, black_list):
            files.append(f)

    return files
def check_file(file, key_words, black_list):
    """
    Check the file against a black list as well as check if it has keywords in
    it.

    Parameters
    ----------
    : file : string
        filename to check
    : key_words : list
        list of strings required to be in file name
    : black_list : list
        list of files that are black listed, that file can't be called
    Returns
    --------
    : : bool
        Returns True if file is not in black_list and contains all key_words
    """
    #check if ile in black list
    if file in black_list:
        return False
    else:
        #check if keyword not in file
        for k in key_words:
            if k not in file:
                return False
        #return true if passes through
        return True

def check_rays(ray_dir, n_rays, fields):
    """
    Check if a directory already contains a given number of trident lrays and
    contains necessary fields

    Parameters
    ----------
    : ray_dir : str
        The path to the directory where rays are held

    : n_rays : int
        The number of lrays that should be in the directory

    : fields : list, str
        List of the fields needed in each lightr ray

    Returns
    --------
    : ray_bool : bool
        `True` if there are `n_rays` in the `ray_dir` and each one contains
        necessary fields. Otherwise returns False
    """

    ray_files = collect_files(ray_dir, key_words=['ray'])

    #check if correct number
    if len(ray_files) == n_rays:
        # check if fields are in each ray

        for rfile in ray_files:
            #load ray file
            try:
                ray = yt.load(f"{ray_dir}/{rfile}")
            except yt.utilities.exceptions.YTOutputNotIdentified:
                print(f"Couldn't load {rfile}. Reconstructing rays")
                return False

            # check each field is in ray
            for fld in fields:
                if ('all', fld) in ray.field_list:
                    pass
                else:
                    print(f"Not all fields present in {rfile}. Reconstructing rays")
                    return False

        # all rays passed
        return True
    else:
        print(f"Only found {len(ray_files)} instead of {n_rays}. Reconstructing rays")
        return False

def combine_astropy_files(directory, kw='ice', outfile=None):

    #get files
    files = collect_files(directory, key_words=['ray', kw])

    tables = []
    # open up tables
    for f in files:
        tables.append(QTable.read(f"{directory}/{f}"))

    if len(tables) >0:
        #combine tables
        main_table = vstack(tables)

        #write table
        if outfile is not None:
            main_table.write(outfile, overwrite=True)
    else:
        out_err = outfile.split('.')[0] + ".out"
        #write out dummy
        f= open(out_err, 'w')
        f.write(f"No files found in {directory} using key_words= ['ray', {kw}]")
        f.close()
        main_table = None
    return main_table

def combine_pandas_files(directory, kw='ice', outfile=None):

    #get files
    files = collect_files(directory, key_words=['ray', kw])

    dfs = []
    # open up tables
    for f in files:
        dfs.append(pd.read_hdf(f"{directory}/{f}"))

    if len(tables) >0:
        #combine tables
        main_table = pd.concat(dfs, ignore_index=True)

        #write table
        if outfile is not None:
            main_table.write_hdf(outfile, mode='w')
    else:
        out_err = outfile.split('.')[0] + ".out"
        #write out dummy
        f= open(out_err, 'w')
        f.write(f"No files found in {directory} using key_words= ['ray', {kw}]")
        f.close()
        main_table = None
    return main_table
