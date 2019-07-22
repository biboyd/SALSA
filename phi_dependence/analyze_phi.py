import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def main():
    #load in data
    ds = argv[1]
    ion=argv[2]
    ion_name = "_".join(ion.split())
    b = int(argv[3])
    scratch = "/mnt/gs18/scratch/users/boydbre1"
    out_dir = "/".join((scratch, f"polar_angle/{ds}/impact_param{b}/{ion_name}"))
    fname = "phi_col_dense_data.npy"
    dat_fname = "/".join((out_dir, fname))
    
    dat=np.load(dat_fname)
    pd_dict = {'Angle':dat[:, 0], 
               'lcd': dat[:,1],
               'azimuth': dat[:, 2]}
    
    df = pd.DataFrame(pd_dict)
    
    means = df.groupby('Angle')['lcd'].mean()
    medians = df.groupby('Angle')['lcd'].median()
    maxes = df.groupby('Angle')['lcd'].max()
    mins = df.groupby('Angle')['lcd'].min()
    
    stats = [maxes, means, medians, mins]
    stat_names = ["max", "mean", "median", "min"]
    fig, ax = plt.subplots()
    for stat, nm in zip(stats, stat_names):
        stat.plot(label=nm)
    ax.set_ylabel(f"Log {ion} Column Density")
    ax.legend()
    
    
    ax.set_title(f'{ion} Col density for b={b} kpc')
    p = np.polyfit(means.index, means.values, 1)
    
    ang = np.linspace(0, 90, 100)
    cd = p[1] + p[0]*ang
    ax.plot(ang, cd, 'k--', label='mean fit')
    
    fig.savefig(f"{out_dir}/stat_plot.png")
    
    #create histogram
    
    az_group = df.groupby('azimuth')
    slopes= az_group.apply(fit_line)
    
    hist_fig, hist_ax = plt.subplots()
    slopes.hist(ax=hist_ax)
    hist_ax.set_title("Slope of lines Fitted to Angle vs Col Dense")
    hist_ax.set_xlabel("Slope")
    hist_ax.set_ylabel('Number')
    
    hist_fig.savefig(f"{out_dir}/slope_hist.png")

def fit_line(df, col1='Angle', col2='lcd', deg=1):
    p = np.polyfit(df[col1].values, df[col2].values, deg)
    
    return p[0]

if __name__ == '__main__':
   main()

