from os.path import join, dirname

import yt
import salsa
import numpy as np
import pandas as pd


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

def test_enzo_generate_catalog(tmp_path, request):

    np.random.seed(567)

    ds = yt.load_sample("IsolatedGalaxy")

    raydir = join(tmp_path, "test_enzo_rays")

    df = salsa.generate_catalog(ds, 2, raydir, ['H I', 'C IV'], method='spice',
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
    key_df = pd.read_csv(join(dirname(request.path), "enzo_key_df.csv"),
                         index_col=0)

    same_catalog(df, key_df)

def test_fire_generate_catalog(tmp_path, request):
    print("tmp_path:", tmp_path)
    np.random.seed(567)

    ds = yt.load_sample("FIRE_M12i_ref11")

    raydir = join(tmp_path, "test_fire_rays")

    c = [29286.1032486 , 31049.29447174, 32589.58339691]
    df = salsa.generate_catalog(ds, 2, raydir, ['H I', 'C IV'], method='spice',
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
    key_df = pd.read_csv(join(dirname(request.path), "fire_key_df.csv"),
                         index_col=0)

    same_catalog(df, key_df)
