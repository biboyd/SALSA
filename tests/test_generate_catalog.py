import yt
import salsa
import numpy as np
import pandas as pd

def test_gen_catalog():

    np.random.seed(567)

    ds = yt.load(f"tests/IsolatedGalaxy/galaxy0030/galaxy0030")
    raydir=f"tests/test_rays/"
    df=salsa.generate_catalog(ds, 2, raydir, ['H I', 'C IV'], method='ice',
                              center=[0.5, 0.5, 0.5], impact_param_lims=(0, 50),
                              ray_length=200,
                              fields=['temperature'],
                              cut_region_filters=["obj['temperature'] < 5e5"],
                              extractor_kwargs={'velocity_res':25, 'absorber_min':12.8, 'frac':0.8},
                              units_dict={'temperature':'K'})

    # check fields were made appropriately
    fields = ['temperature', 'H_p0_number_density', 'C_p3_number_density']
    assert salsa.utils.check_rays(raydir, 2, fields)

    # load in dataframe to compare to
    key_df = pd.read_csv(f"tests/key_df.csv", index_col=0)

    # check number of columns/rows
    assert len(df) == len(key_df)
    assert len(df.columns) == len(key_df.columns)

    # test extraction values same as key
    for i in range(len(df)):
        # extract correct row to compare to
        ion=df.loc[i, 'name']
        index = df.loc[i, 'absorber_index']
        compare = key_df.query(f'(name == "{ion}") & (absorber_index == "{index}")').iloc[0]

        # compare column densities
        cd_diff = abs(df.loc[i, 'col_dens'] - compare['col_dens'])/compare['col_dens']
        assert cd_diff < 1e-4

        #compare temperatures
        temp_diff = abs(df.loc[i, 'temperature'] - compare['temperature'])/compare['temperature']
        assert temp_diff < 1e-4
