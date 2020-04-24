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

def double_hist(data1, data2, bins, hist_range=None, ax=None, color1='tab:blue', color2='tab:orange',label1='Spectacle', label2="rCloudy"):
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

## MAKE COLUMN DENSITY HISTO comparing Spectacle and rCloudy
def col_dense_histo(spect_df, cloudy_df, z,ion="O VI", bins=20, hist_range=None, outdir='./'):
    ion_u=f"{ion.split()[0]}_{ion.split()[1]}"
    spect_cd = spect_df['col density']
    cloudy_cd = cloudy_df['col density']

    fig, ax = plt.subplots(1)
    double_hist(spect_cd, cloudy_cd, bins, hist_range, ax=ax, color1='tab:purple', color2='black') 

    ax.set_title(f"{ion} Column Density for Spectacle and rCloudy (z={z:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel(f"Log( $N_{ {ion} }$ )")
    ax.legend()

    fig.savefig(f"{outdir}/col_density_{ion_u}_{z:0.2f}.png", dpi=400)

## MAKE VELOCITY HISTO
def los_velocity_histo(spect_df, cloudy_df, z,ion="O VI", bins=25, hist_range=None,outdir='./'):
    ion_u=f"{ion.split()[0]}_{ion.split()[1]}"
    spect_vel = spect_df['avg_velocity']
    cloudy_vel = cloudy_df['avg_velocity']

    fig, ax = plt.subplots(1)
    double_hist(spect_vel, cloudy_vel, bins, hist_range, ax=ax, color1='tab:purple', color2='black') 

    ax.set_title(f"{ion} Velocity for Spectacle and rCloudy (z={z:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel("Line of Sight Velocity (km/s)")
    ax.legend()

    fig.savefig(f"{outdir}/velocity_{ion_u}_{z:0.2f}.png", dpi=400)

#### make col density vs b scatter plot

#MAKE b histogram
def b_histo(spect_df, z, b_lim=None, bins=15, ion="O VI",outdir='./', filtered=False):
    ion_u=f"{ion.split()[0]}_{ion.split()[1]}"
    fld= 'vel_doppler'
    if filtered is True:
        #extract inflow and outflow information
        inflow = spect_df[ spect_df['flow'] == 'Inflow'][fld]
        outflow = spect_df[ spect_df['flow'] == 'Outflow'][fld]
        fig, ax = plt.subplots(1)

        double_hist(inflow, outflow, bins, b_lim, ax=ax, label1='Inflow', label2='Outflow')

    else:
        data = spect_df[fld]
        fig, ax = plt.subplots(1)
        data_histo, data_edges= np.histogram(data, range=b_lim, bins=bins)
        data_histo = data_histo/data.size *100
        ax.hist(data_edges[:-1], data_edges, weights=data_histo, label="# absorbers: {data.size}")
    
    ax.set_title(f"{ion} Doppler Velocity (z={z:0.2f})")
    ax.set_xlabel("b (km/s)")
    ax.set_ylabel("% of absorbers")
    ax.set_ylim(0, 50)
    ax.legend()
    fig.savefig(f"{outdir}/b_inflow_outflow_{ion_u}_{z:0.2f}.png", dpi=400)

## MAKE VELOCITY HISTO
def inflow_outflow_histos(cloudy_df, z,ion="O VI", bins=15, hist_range=None, outdir='./'):
    ion_u=f"{ion.split()[0]}_{ion.split()[1]}"

    titles=['Metallicity', 'Radius', 'Density', 'Temperature']
    fields=['log_avg_metallicity', 'radius', 'log_avg_density', 'log_avg_temperature']
    xlabels=["Log( $Z_{\odot}$ )", 'Radius (kpc)', 'Log Density (g/cm^3)', 'Log Temperature (K)'] 

    for title, fld, xlabel, curr_range in zip(titles, fields, xlabels, hist_range):
        #extract inflow metal and outflow metal information
        inflow = cloudy_df[ cloudy_df['flow'] == 'Inflow'][fld]
        outflow = cloudy_df[ cloudy_df['flow'] == 'Outflow'][fld]

        fig, ax = plt.subplots(1)
        double_hist(inflow, outflow, bins, curr_range, ax=ax, label1='Inflow', label2='Outflow')

        ax.set_title(f"{ion} {title} for Inflows and Outflows (z={z:0.2f})")
        ax.set_ylabel("% of absorbers")
        ax.set_xlabel(xlabel)
        ax.set_ylim(0, 50)
        ax.legend(loc='upper left')

        fig.savefig(f"{outdir}/{title}_inflow_outflow_{ion_u}_{z:0.2f}.png", dpi=400)

# run plot functions

if __name__ == '__main__':
    
    #define data set info
    dsname=sys.argv[1]
    max_b=200
    frac=0.8
    ion=sys.argv[2]
    

    ion_u = "_".join(ion.split(" "))
    cgm='CGM'
    ds_fname = f"/mnt/gs18/scratch/users/boydbre1/cosmological/foggie/{dsname}/{dsname}"
    if len(sys.argv) == 3:
        arr_dir=f"/mnt/gs18/scratch/users/boydbre1/analysis/foggie/{dsname}/max_impact{max_b}/output_{cgm}{frac}/ion_{ion_u}"
        ray_dir=f"/mnt/gs18/scratch/users/boydbre1/analysis/foggie/{dsname}/max_impact{max_b}/rays"

        #define the outdirectory
        outdir=f"plots/{ion_u}_plots"
        makedirs(outdir, exist_ok=True)
        filtered=False

    elif len(sys.argv) ==4 and sys.argv[3] == 'filtered':
        arr_dir=f"/mnt/home/boydbre1/data/O_VI_analysis/{dsname}_flow/ion_{ion_u}"

        ray_dir=f"/mnt/gs18/scratch/users/boydbre1/analysis/foggie/{dsname}/max_impact{max_b}/rays"

        #define the outdirectory
        outdir=f"plots_filtered/{ion_u}_plots"
        makedirs(outdir, exist_ok=True)
        filtered=True

    #load in files
    absorber_array=np.load(f"{arr_dir}/all_absorbers.npy")
    ds = yt.load(ds_fname)
    all_impact_param = np.load(f"{ray_dir}/impact_parameter.npy")
    header = np.load(f"{arr_dir}/absorber_info_header.npy")
    df = pd.DataFrame(data=absorber_array, columns=header)
    z = ds.current_redshift

    #divide data into spectacle and cloudy
    spect_df = df[ df['spectacle'] == 1.]
    cloudy_df = df[ df['cloudy'] == 1.]
    print('z:',z, 
          '\nspec #:',spect_df['ray_num'].count(),
          '\nrcloud #:',cloudy_df['ray_num'].count(),
          )

    if filtered is False:
        # DEFINE INFLOW OUTLFOW STUFF
        cloudy_df['flow'] = cloudy_df['radial_velocity'].apply(lambda r : "Inflow" if (r < 0) else "Outflow" )
    elif filtered is True:
        #set flow for already filtered data
        cloudy_df['flow'] = cloudy_df['inflow'].apply(lambda i : "Inflow" if (i == 1.) else "Outflow")
        spect_df['flow'] = spect_df['inflow'].apply(lambda i : "Inflow" if (i == 1.) else "Outflow")

    # LOG SOME VARS FOR PAIRPLOT
    variables =['avg_density', 'avg_metallicity', 'avg_temperature']
    for v in variables:
        logv = 'log_' + v
        cloudy_df[logv] = np.log10(df[v])
    log_var = ['col density', 'log_avg_density', 'log_avg_temperature', 'log_avg_metallicity', 'radius', 'radial_velocity']

    # labels to use in final plot instead
    axis_labels={'log_avg_density':"Log( Density ) ($g/cm^3$)", 
                      'log_avg_metallicity':"Log( Metallicity ) ($Z_{\odot}$)", 
                      'radius':'Radial Distance ($kpc$)', 
                      'log_avg_temperature':"Log( Temperature ) ($K$)", 
                      'radial_velocity': "Radial Velocity ($km/s$)"}

    #create pairplot
    # MAKE AND SAVE PAIRPLOT
    sns.set_palette('colorblind')
    pp= sns.pairplot(cloudy_df, vars=log_var, hue='flow', markers='o', diag_kind='hist', diag_kws=dict(alpha=0.7), plot_kws=dict(alpha=0.7))
    pp.fig.suptitle(f"{ion} properties for z= {z:0.2f}",size=18, y=1.05)
   
    # change the axis labels
    n_var=len(log_var)
    for i in range(n_var):
        for j in range(n_var):
            #find x/y label
            xlabel = pp.axes[i][j].get_xlabel()
            ylabel = pp.axes[i][j].get_ylabel()

            # replace if found in labels dict
            if xlabel in axis_labels.keys():
                pp.axes[i][j].set_xlabel(axis_labels[xlabel])

            if ylabel in axis_labels.keys():
                pp.axes[i][j].set_ylabel(axis_labels[ylabel])
            
    # save pairplot 
    pp.savefig(f"{outdir}/pairplot_{ion_u}_{z:0.2f}.png")
    
    
    hist_range_dict = {"O VI":[(13., 16.), (-300, 200), (0, 150), (-2, 0.05), (-28.75, -25), (4, 7)],
                       "H II":None
                      }
    if ion in hist_range_dict.keys():
        cd_lim, los_vel_lim, b_lim, Z_lim, dense_lim, temp_lim = hist_range_dict[ion]
    else:
        cd_lim = None
        los_vel_lim=None
        b_lim = None
        Z_lim = None
    #run plot functions
    print(f"creating plotting stuff rn, goin into {outdir}")
    col_dense_histo(spect_df, cloudy_df, z, ion=ion,  hist_range=cd_lim, outdir=outdir)
    los_velocity_histo(spect_df, cloudy_df, z, ion=ion,  hist_range=los_vel_lim, outdir=outdir)
    b_histo(spect_df, z, ion=ion, b_lim=b_lim,outdir=outdir, filtered=filtered)
    inflow_outflow_histos(cloudy_df, z, ion=ion, bins=10,hist_range=(Z_lim, (0, 200), dense_lim, temp_lim),outdir=outdir)
    
