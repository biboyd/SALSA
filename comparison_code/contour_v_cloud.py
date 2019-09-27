# Brendan Boyd 09/26/19
# Compare two methods of extracting features from number density
# contour method by Hillary Egan and cloud method by Jason Tumlinson
#

import matplotlib.pyplot as plt
import numpy as np
import yt
import trident
from matplotlib.patches import Patch
from os import makedirs
import sys
sys.path.insert(0, "/mnt/home/boydbre1/Repo/CGM/multi_plot_movie")
sys.path.insert(0, "/home/bb/Repo/CGM/multi_plot_movie")

from multi_plot import multi_plot
from center_finder import find_center


def main_compare(ds_filename, ray_filename, ion_list, out_dir="./"):
    """
    plot the intervals of the two methods of extraction on top of
    a number density plot vs ray length

    Parameters:
        ds_filename :string: path to dataset
        ray_filename :string: path to the light ray file
        ion_list :list: List of ion names to look at (ie "H I" "C IV"...)
        out_dir :string: directory in which to save the plots
    Returns:
        none
    """
    makedirs(out_dir, exist_ok=True)
    center, nvec, redshift, bv = find_center(ds_filename)
    for ion in ion_list:
        mp = multi_plot(ds_filename, ray_filename,
                        ion_name=ion,
                        redshift=redshift)

        #create plot of number density
        fig, ax = plt.subplots()
        mp.plot_num_dense_los_vel(ax_num_dense=ax)

        #extract features using two methods
        ctour_intervals, ctour_lcd = mp.get_contour_intervals()
        cloud_intervals, cloud_lcd = mp.get_cloud_intervals()
        l_list = mp.ray.data['l'].in_units('kpc')

        #plot methods
        ctour_args=dict(color='tab:blue', alpha=0.5)
        cloud_args=dict(color='tab:orange', alpha=0.5)

        ctour_tot= plot_intervals(ax, l_list, ctour_intervals, ctour_lcd, plot_args=ctour_args)
        cloud_tot = plot_intervals(ax, l_list, cloud_intervals, cloud_lcd, plot_args=cloud_args)

        ctour_patch = Patch(label=f"contour: {ctour_tot:.1f}", **ctour_args)
        cloud_patch = Patch(label=f"cloud: {cloud_tot:.1f}", **cloud_args)

        ax.legend(handles=[ctour_patch, cloud_patch])

        outfile=f"{mp.ion_p_name():}_feature_comp.png"
        fig.savefig(f"{out_dir}/{outfile}")

def plot_intervals(ax, l_list, intervals, lcd_list, plot_args=dict(color='tab:grey', alpha=0.5)):
    tot_lcd=0
    for i in range(len(intervals)):
        b, e = intervals[i]
        curr_lcd = lcd_list[i]

        #plot interval
        ax.axvspan(l_list[b], l_list[e], **plot_args)
        tot_lcd += 10**curr_lcd
    return np.log10(tot_lcd)

if __name__ == '__main__':
    ds =sys.argv[1]
    ray = sys.argv[2]
    ion = ['H I']
    main_compare(ds, ray, ion)
