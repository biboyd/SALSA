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


def load_files(ds, refinement='cool', ion="O VI", cut1="cgm", cut2=None):
    ds_fname = f"/mnt/gs18/scratch/users/boydbre1/cosmological/{refinement}_refinement/{ds}/{ds}"

    ion_u = "_".join(ion.split(" "))
    #load data
    dataDir= f"{homeDir}/data/absorber_data/{refinement}_refinement/ion_{ion_u}/"

    # read in cut1 file
    cut1_u = "_".join(cut1.split(" "))
    cut1_file = f"{dataDir}/{cut1_u}/{ds}_absorbers.h5"
    cut1_table = Table.read(cut1_file)

    cut1_table['cuts'] = cut_alias_dict[cut1_u]

    #read in cut1 file
    if cut2 is None:
        final_table = cut1_table

    else:
        cut2_u = "_".join(cut2.split(" "))
        cut2_file = f"{dataDir}/{cut2_u}/{ds}_absorbers.h5"
        cut2_table = Table.read(cut2_file)
        cut2_table['cuts'] = cut_alias_dict[cut2_u]

        final_table = vstack([cut1_table, cut2_table])

    #try to convert to pandas
    try:
        df = final_table.to_pandas()
        return df
    except ImportError:
        print("Pandas not installed. returning py table instead")
        return final_table
