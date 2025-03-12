import os

import yt
import salsa
import numpy as np
import pandas as pd

def test_enzo_generate_catalog():

    np.random.seed(567)

    ds = yt.load_sample("IsolatedGalaxy")
    raydir=f"tests/test_enzo_rays/"

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
    key_df = pd.read_csv(f"tests/enzo_key_df.csv", index_col=0)

    same_catalog(df, key_df)

def test_fire_generate_catalog():

    np.random.seed(567)

    ds = yt.load_sample("FIRE_M12i_ref11")
    raydir=f"tests/test_fire_rays/"

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
    key_df = pd.read_csv(f"tests/fire_key_df.csv", index_col=0)

    df.to_csv('tests/fire_df.csv')
    same_catalog(df, key_df)

def clear_ray_dir(folder):
    """
    Explicitly remove files to make sure weird extra stuff isn't being added.
    """
    for i in range(2):
        try:
            os.remove(os.path.join(folder, f"ray{i}.h5"))
        except FileNotFoundError:
            continue
    try:
        os.remove("impact_parameters.npy")
    except FileNotFoundError:
        pass

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
