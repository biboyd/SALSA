import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt
import sys
sys.path.insert(0, "/mnt/home/boydbre1/Repo/CGM/multi_plot_movie/")
import multi_plot
from center_finder import find_center

from os import makedirs
import seaborn as sns


# plot functions

def double_hist(data1, data2, bins, hist_range=None, ax=None, label1='Spectacle', label2="rCloudy"):
    if ax is None:
        fig, ax = plt.subplots(1)
    else:
        fig = ax.figure
    
    data1_histo, data1_edges= np.histogram(data1, range=hist_range, bins=bins)
    data1_histo = data1_histo/data1.size *100

    ax.hist(data1_edges[:-1], data1_edges, weights=data1_histo, label=f"{label1}: {data1.size}")

    data2_histo, data2_edges= np.histogram(data2, range=hist_range, bins=bins)
    data2_histo = data2_histo/data2.size *100

    ax.hist(data2_edges[:-1], data2_edges, weights=data2_histo, histtype='step', label=f"{label2}: {data2.size}")
    
    return fig, ax

## MAKE COLUMN DENSITY HISTO
def col_dense_histo(spect_df, cloudy_df, ds, outdir='./'):
    bins=20
    hist_range=(13, 14.75)
    spect_cd = spect_df['col density']
    cloudy_cd = cloudy_df['col density']

    fig, ax = plt.subplots(1)
    double_hist(spect_cd, cloudy_cd, bins, hist_range, ax=ax)

    ax.set_title(f"O VI Column Density for Spectacle and rCloudy (z={ds.current_redshift:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel("Log( $N_{OVI}$ )")
    ax.legend()

    fig.savefig(f"{outdir}/col_density_O_VI_{ds.current_redshift:0.2f}.png", dpi=400)

## MAKE VELOCITY HISTO
def los_velocity_histo(spect_df, cloudy_df, ds, outdir='./'):
    bins=25
    hist_range=(-200, 200)
    spect_vel = spect_df['avg_velocity']
    cloudy_vel = cloudy_df['avg_velocity']

    fig, ax = plt.subplots(1)
    double_hist(spect_vel, cloudy_vel, bins, hist_range, ax=ax)

    ax.set_title(f"O VI Velocity for Spectacle and rCloudy (z={ds.current_redshift:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel("Line of Sight Velocity (km/s)")
    ax.legend()

    fig.savefig(f"{outdir}/velocity_O_VI_{ds.current_redshift:0.2f}.png", dpi=400)

#### make col density vs b scatter plot

#MAKE COL DENSITY VS b SCATTER
def cd_vs_b_histo(spect_df, cloudy_df, ds, outdir='./'):
    spect_cd = spect_df['col density']
    spect_b = spect_df['vel_doppler']

    fig, ax = plt.subplots(1)

    ax.scatter(spect_b, spect_cd, marker='.')

    ax.set_title(f"O VI Column Density vs b (z={ds.current_redshift:0.2f})")
    ax.set_xlabel("b (km/s)")
    ax.set_ylabel("Log( $N_{OVI}$ )")
    fig.savefig(f"{outdir}/cd_vs_b_O_VI_{ds.current_redshift:0.2f}.png", dpi=400)

## MAKE VELOCITY HISTO
def inflow_outflow_histo(cloudy_df, ds, outdir='./'):

    bins=15
    hist_range=None#(-200, 200)

    #extract inflow metal and outflow metal information
    inflow_metal = cloudy_df[ cloudy_df['flow'] == 'Inflow'].log_avg_metallicity
    outflow_metal = cloudy_df[ cloudy_df['flow'] == 'Outflow'].log_avg_metallicity

    fig, ax = plt.subplots(1)
    double_hist(inflow_metal, outflow_metal, bins, hist_range, ax=ax, label1='Inflow', label2='Outflow')

    ax.set_title(f"O VI Metallicity for Inflows and Outflows (z={ds.current_redshift:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel("Log( $Z_\odot$ )")
    ax.legend(loc='upper left')

    fig.savefig(f"{outdir}/metallcity_inflow_outflow_O_VI_{ds.current_redshift:0.2f}.png", dpi=400)


# run plot functions

if __name__ == '__main__':
    
    #define data set info
    dsname=sys.argv[1]
    max_b=200
    frac=0.8
    ion="O_VI"
    cgm='CGM'
    ds_fname = f"/mnt/gs18/scratch/users/boydbre1/cosmological/foggie/{dsname}/{dsname}"
    arr_dir=f"/mnt/gs18/scratch/users/boydbre1/metallicity_analysis/foggie/{dsname}/max_impact{max_b}/output_{cgm}{frac}/ion_{ion}"
    ray_dir=f"/mnt/gs18/scratch/users/boydbre1/metallicity_analysis/foggie/{dsname}/max_impact{max_b}/rays"


    #load in files
    absorber_array=np.load(f"{arr_dir}/all_absorbers.npy")
    ds = yt.load(ds_fname)
    all_impact_param = np.load(f"{ray_dir}/impact_parameter.npy")
    header = np.load(f"{arr_dir}/absorber_info_header.npy")
    df = pd.DataFrame(data=absorber_array, columns=header)

    #divide data into spectacle and cloudy
    spect_df = df[ df['spectacle'] == 1.]
    cloudy_df = df[ df['cloudy'] == 1.]
    print('z:',ds.current_redshift, 
          '\nspec #:',spect_df['ray_num'].count(),
          '\nrcloud #:',cloudy_df['ray_num'].count(),
          )

    # DEFINE INFLOW OUTLFOW STUFF
    cloudy_df['flow'] = cloudy_df['radial_velocity'].apply(lambda r : "Inflow" if (r < 0) else "Outflow" )
    # LOG SOME VARS FOR PAIRPLOT
    variables =['avg_density', 'avg_metallicity', 'avg_temperature']
    for v in variables:
        logv = 'log_' + v
        cloudy_df[logv] = np.log10(df[v])
    log_var = ['col density','log_avg_density', 'log_avg_metallicity', 'avg_velocity', 'log_avg_temperature', 'radial_velocity']

    #define the outdirectory
    outdir="ovi_plots"
    makedirs(outdir, exist_ok=True)

    #create pairplot
    # MAKE AND SAVE PAIRPLOT
    pp= sns.pairplot(cloudy_df, vars=log_var, hue='flow', markers='o')
    pp.savefig(f"{outdir}/pairplot_O_VI_{ds.current_redshift:0.2f}.png", dpi=400)
    
    #run plot functions
    print(f"creating plotting stuff rn, goin into {outdir}")
    col_dense_histo(spect_df, cloudy_df, ds, outdir=outdir)
    los_velocity_histo(spect_df, cloudy_df, ds, outdir=outdir)
    cd_vs_b_histo(spect_df, cloudy_df, ds, outdir=outdir)
    inflow_outflow_histo(cloudy_df, ds, outdir=outdir)
