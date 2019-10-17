# Brendan Boyd 09/26/19
# Compare two methods of extracting features from number density
# contour method by Hillary Egan and cloud method by Jason Tumlinson
#
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt
import trident
from matplotlib.patches import Patch
from functools import reduce
from os import makedirs
import sys
sys.path.insert(0, "/mnt/home/boydbre1/Repo/CGM/multi_plot_movie")
sys.path.insert(0, "/home/bb/Repo/CGM/multi_plot_movie")

from multi_plot import multi_plot
from center_finder import find_center


def plot_compare(mp, out_dir="./", fname=None, char_frac=0.5, dense_frac=0.85):
    """
    plot the intervals of the two methods of extraction on top of
    a number density plot vs ray length

    Parameters:
        ds_fname :string: path to dataset
        ray_fname :string: path to the light ray file
        ion_list :list: List of ion names to look at (ie "H I" "C IV"...)
        out_dir :string: directory in which to save the plots
    Returns:
        none
    """
    makedirs(out_dir, exist_ok=True)
    for ion in mp.ion_list:

        #create plot of number density
        fig, ax = plt.subplots()
        mp.plot_num_dense_los_vel(ax_num_dense=ax)

        #extract features using two methods
        ctour_intervals, ctour_lcd = mp.get_contour_intervals(char_density_frac=char_frac, force=True)
        cloud_intervals, cloud_lcd = mp.get_cloud_intervals(coldens_fraction=dense_frac)
        l_list = mp.ray.data['l'].in_units('kpc')

        #plot methods
        ctour_args=dict(color='tab:blue', alpha=0.5)
        cloud_args=dict(color='tab:orange', alpha=0.5)

        ctour_tot= plot_intervals(ax, l_list, ctour_intervals, ctour_lcd, plot_args=ctour_args)
        cloud_tot = plot_intervals(ax, l_list, cloud_intervals, cloud_lcd, plot_args=cloud_args)

        ctour_patch = Patch(label=f"contour: {ctour_tot:.1f}", **ctour_args)
        cloud_patch = Patch(label=f"cloud: {cloud_tot:.1f}", **cloud_args)

        ax.legend(handles=[ctour_patch, cloud_patch])

        if fname is None:
            outfile=f"{mp.ion_p_name():}_feature_comp.png"
        else:
            outfile = fname
        fig.savefig(f"{out_dir}/{outfile}")
        plt.close(fig)

def plot_intervals(ax, l_list, intervals, lcd_list, plot_args=dict(color='tab:grey', alpha=0.5)):
    tot_lcd=0
    for i in range(len(intervals)):
        b, e = intervals[i]
        curr_lcd = lcd_list[i]

        #plot interval
        ax.axvspan(l_list[b], l_list[e], **plot_args)
        tot_lcd += 10**curr_lcd
    return np.log10(tot_lcd)

def get_info(mp, char_frac=0.5, dense_frac=0.85):

    ctour_intervals, ctour_lcd = mp.get_contour_intervals(char_density_frac=char_frac, force=True)
    cloud_intervals, cloud_lcd = mp.get_cloud_intervals(coldens_fraction=dense_frac)

    tot_ctour = 0
    for lcd in ctour_lcd: 
        tot_ctour += 10**lcd
    tot_ctour_lcd = np.log10(tot_ctour)
    tot_cloud = 0
    for lcd in cloud_lcd: 
        tot_cloud += 10**lcd
    tot_cloud_lcd = np.log10(tot_ctour)
    
    n_contour = len(ctour_intervals)
    n_cloud = len(cloud_intervals)

    """
    ADD A SIMILAR THING THAT CHECKS FOR PERCENT OF NUMBER DENSITY THEY OVERLAP
    """
    #check for overlaps
    length_arr = np.array(mp.ray.data["l"])
    ctour_arrs=np.array([])
    cloud_arrs=np.array([])
    for i, f in ctour_intervals:
        ctour_arrs = np.concatenate((ctour_arrs,length_arr[i:f]), 0)
    for i, f in cloud_intervals:
        cloud_arrs = np.concatenate((cloud_arrs,length_arr[i:f]), 0)

    intersect_arr = np.intersect1d(ctour_arrs, cloud_arrs)
    num_intersect = len(intersect_arr)
    num_covered = len(ctour_arrs) + len(cloud_arrs)
    if num_covered > 0:
        percent_intersect = num_intersect/num_covered * 100
    else:
        percent_intersect = 0

    print(f"percent intersect {percent_intersect: .2f}%")
    return (n_contour, n_cloud, tot_ctour_lcd, tot_cloud_lcd, percent_intersect)


def parameter_sweep(mp, steps=25, out_dir='./', char_bound=[0.1, 0.9],dense_bound=[0.5, 0.95]):
    dense_frac_array=np.linspace(dense_bound[0], dense_bound[1], steps)
    char_frac_array= np.linspace(char_bound[0], char_bound[1], steps)

    all_info_array = np.empty((25**2, 7))
    for i in range(steps):
        char_frac = char_frac_array[i]
        for j in range(steps):
            dense_frac = dense_frac_array[j]

            #create plots
            out_fname = f"char_{char_frac:0.2f}_dense_{dense_frac:0.2f}_comp.png"
            plot_compare(mp, out_dir=out_dir, fname=out_fname, char_frac=char_frac, dense_frac=dense_frac)
            #get info on overlap
            curr_info = get_info(mp, char_frac=char_frac, dense_frac=dense_frac)
            #fill info into array
            all_info_array[25*i + j, 0] = char_frac
            all_info_array[25*i + j, 1] = dense_frac
            #print(curr_info)
            #print(all_info_array[25*i + j, 2:])
            all_info_array[25*i + j, 2:] = curr_info

    np.save(f"{out_dir}/all_info.npy", all_info_array)

if __name__ == '__main__':
    ds =sys.argv[1]
    ray = sys.argv[2]
    ion = sys.argv[3]
    out_dir = sys.argv[4]
    #remove path and ".h5" from ray
    ray_name = ray.split("/")[-1]
    ray_name = ray_name.split('.')[0]
    out_dir = f"{out_dir}/{ray_name}_{ion[0]}_{ion[-1]}"
    c, n, r, bv = find_center(ds)
    mp = multi_plot(ds, ray, ion_name=ion, redshift=r)
    parameter_sweep(mp, steps=25, out_dir=out_dir)
