import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt
import sys
from astropy.table import Table, vstack
from os import makedirs
import seaborn as sns

homeDir="/mnt/home/boydbre1"

from CGM.general_utils.filter_definitions import cut_alias_dict
from CGM.general_utils.center_finder import find_center

def double_hist(data1, data2, bins, hist_range=None, ax=None, color1='tab:blue', color2='tab:orange',label1='Spectacle', label2="ICE"):
    if ax is None:
        fig, ax = plt.subplots(1)
    else:
        fig = ax.figure

    data1_histo, data1_edges= np.histogram(data1, range=hist_range, bins=bins)
    data1_histo = data1_histo/data1.size *100

    ax.hist(data1_edges[:-1], data1_edges, weights=data1_histo, color=color1, label=f"{label1}: {data1.size}")

    data2_histo, data2_edges= np.histogram(data2, range=hist_range, bins=bins)
    data2_histo = data2_histo/data2.size *100

    ax.hist(data2_edges[:-1], data2_edges, weights=data2_histo, histtype='step',color=color2, label=f"{label2}: {data2.size}")

    return fig, ax


def load_files(ds, refinement='cool', ion="O VI", cuts=["cgm"]):
    ds_fname = f"/mnt/gs18/scratch/users/boydbre1/cosmological/{refinement}_refinement/{ds}/{ds}"

    ion_u = "_".join(ion.split(" "))
    #load data
    dataDir= f"{homeDir}/data/absorber_data/{refinement}_refinement/ion_{ion_u}/"

    # load tables
    table_list=[]
    for c in cuts:
        c_u = "_".join(c.split(" "))
        c_file = f"{dataDir}/{c_u}/{ds}_absorbers.h5"
        try: 
            c_table = Table.read(c_file)
            c_table['cuts'] = cut_alias_dict[c_u]

            table_list.append(c_table)
        except OSError:
            print(f"couldn't load {c_file}")

    if len(table_list) == 0:
        raise FileNotFoundError("Could not find any h5 files to load")
    else:
        #combine tables 
        final_table = vstack(table_list)

    #try to convert to pandas
    try:
        df = final_table.to_pandas()
        return df
    except ImportError:
        print("Pandas not installed. returning py table instead")
        return final_table
