import os
import sys

import yt
import salsa
import numpy as np
import pandas as pd

def test_enzo_generate_catalog(testdir):

    np.random.seed(567)

    ds = yt.load_sample("IsolatedGalaxy")

    raydir=os.path.join(testdir, "test_enzo_rays")
    clear_ray_dir(raydir)

    df=salsa.generate_catalog(ds, 2, raydir, ['H I', 'C IV'], method='spice',
                              center=[0.5, 0.5, 0.5], impact_param_lims=(0, 50),
                              ray_length=200,
                              fields=['temperature'],
                              ftype="gas",
                              extractor_kwargs={'velocity_res':25, 'absorber_min':12.8, 'frac':0.8},
                              units_dict={'temperature':'K'})

    # check fields were made appropriately
    fields = ['temperature', 'H_p0_number_density', 'C_p3_number_density']
    assert salsa.utils.check_rays(raydir, 2, fields)

    # load in dataframe to compare to
    key_df = pd.read_csv(os.path.join(testdir, "enzo_key_df.csv"), index_col=0)

    same_catalog(df, key_df)

def test_fire_generate_catalog(testdir):

    np.random.seed(567)

    ds = yt.load_sample("FIRE_M12i_ref11")

    raydir=os.path.join(testdir,"test_fire_rays")
    clear_ray_dir(raydir)

    c=[29286.1032486 , 31049.29447174, 32589.58339691]
    df=salsa.generate_catalog(ds, 2, raydir, ['H I', 'C IV'], method='spice',
                              center=c, impact_param_lims=(0, 50),
                              ray_length=200,
                              fields=['temperature'],
                              ftype='PartType0',
                              extractor_kwargs={'velocity_res':25, 'absorber_min':12.8, 'frac':0.8},
                              units_dict={'temperature':'K'})

    # check fields were made appropriately
    fields = ['temperature', 'H_p0_number_density', 'C_p3_number_density']
    assert salsa.utils.check_rays(raydir, 2, fields)

    # load in dataframe to compare to
    key_df = pd.read_csv(os.path.join(testdir, "fire_key_df.csv"), index_col=0)

    same_catalog(df, key_df)

def clear_ray_dir(folder):
    """
    Clear out the ray directory so that rays will be re-generated even if
    previous tests have been run. 
    Explicitly remove each file to make sure weird extra stuff isn't being added,
    and then delete the folder to make sure it is fully empty.
    """
    for i in range(2):
        try:
            os.remove(os.path.join(folder, f"ray{i}.h5"))
        except FileNotFoundError:
            continue # tests may not have been run before
    try:
        os.remove(os.path.join(folder, "impact_parameter.npy"))
    except FileNotFoundError:
        pass # tests may not have been run before

    try:
        os.rmdir(folder)
    except FileNotFoundError: # as opposed to OSError, where dir is not empty
        pass # tests may not have been run before

def same_catalog(df, key_df):
    # check number of columns/rows
    assert len(df) == len(key_df)
    assert len(df.columns) == len(key_df.columns)

    # test extraction values same as key
    for i in range(len(df)):
        # compare column densities
        cd_diff = abs(df.loc[i, 'col_dens'] - key_df.loc[i, 'col_dens'])/key_df.loc[i, 'col_dens']
        assert cd_diff < 1e-4

        #compare temperatures
        temp_diff = abs(df.loc[i, 'temperature'] - key_df.loc[i,'temperature'])/key_df.loc[i,'temperature']
        assert temp_diff < 1e-4

if __name__ == "__main__":
    ray_basedir = os.path.dirname(sys.argv[0])  # path to this executable
    test_enzo_generate_catalog(ray_basedir)
    test_fire_generate_catalog(ray_basedir)
    print("Tests successful.")