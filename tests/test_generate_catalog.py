from os.path import join, dirname

import yt
import salsa
import numpy as np
from astropy.io import ascii

def same_catalog(df, key_df):
    # check number of columns/rows
    assert len(df) == len(key_df)
    assert len(df.columns) == len(key_df.columns)

    # compare column densities
    cd_diff = abs(df['col_dens'] - key_df['col_dens'])/key_df['col_dens']
    assert (cd_diff.value < 1e-4).all()

        #compare temperatures
    temp_diff = abs(df['temperature'] - key_df['temperature']) \
        / key_df['temperature']
    assert (temp_diff.value < 1e-4).all()

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
    key_df = ascii.read(join(dirname(request.path), "enzo_key.ecsv"))

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
    key_df = ascii.read(join(dirname(request.path), "fire_key.ecsv"))

    same_catalog(df, key_df)
