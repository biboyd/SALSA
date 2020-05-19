import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt
import argparse

from os import makedirs
import seaborn as sns

from CGM.absorber_analysis.helper_functions import double_hist, load_files, homeDir
from CGM.general_utils.filter_definitions import cut_alias_dict, axis_labels_dict, hist_range_dict

def main(df, ion, cuts, outdir):
    #load dataframe
    z=df['redshift'][0]

    # LOG SOME VARS FOR PAIRPLOT
    variables =['density', 'metallicity', 'temperature']
    for v in variables:
        logv = 'log_' + v
        df[logv] = np.log10(df[v])

    plot_var = ['col_dens', 'log_density', 'log_temperature', 'log_metallicity', 'radius', 'radial_velocity']

    #plot pairplot
    pp =create_pairplot(df, plot_var=plot_var)
    pp.savefig(f"{outdir}/pairplot_{z:0.2f}.png")

    #make double hist if just two cuts
    if len(cuts)==2:
        # check general_utils for dictionary of limits
        if ion in hist_range_dict:
            my_hist_dict=hist_range_dict[ion]
        else:
            #empty dictionary
            my_hist_dict={}

        #plot histograms of plot_vars
        for var in plot_var:
            if var in my_hist_dict.keys():
                hist_range = my_hist_dict[var]
            else:
                hist_range=None

            fig = plot_histograms(df, var, ion=ion, bins=15, hist_range=hist_range)
            fig.savefig(f"{outdir}/{var}_{z:0.2f}.png", dpi=400)
            plt.close(fig)

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


def plot_histograms(df, var, ion="O VI", bins=15, hist_range=None, outdir='./'):

    #default hist_range to min/max
    if hist_range is None:
        hist_range = [df[var].min(), df[var].max()]

    #return cut titles
    cut1, cut2 = df['cuts'].unique()
    cut1_dat = df[ df['cuts'] == cut1 ][var]
    cut2_dat = df[ df['cuts'] == cut2 ][var]
    z = df['redshift'][0]

    if var in axis_labels_dict.keys():
        xlabel=axis_labels_dict[var]
    else:
        xlabel=var

    fig, ax = plt.subplots(1)
    double_hist(cut1_dat, cut2_dat, bins, hist_range, ax=ax, label1=cut1, label2=cut2)
    ax.set_title(f"{ion} {var} for {cut1} and {cut2} (z={z:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, 50)
    ax.legend(loc='upper left')

    return fig

def create_pairplot(df, plot_var):

    #create pairplot
    # MAKE AND SAVE PAIRPLOT
    sns.set_palette('colorblind')
    pp= sns.pairplot(df, vars=plot_var, hue='cuts', markers='o', diag_kind='hist', corner=True, diag_kws=dict(alpha=0.7), plot_kws=dict(alpha=0.7))
    pp.fig.suptitle(f"{ion} properties for z= {df['redshift'][0]:0.2f}",size=18, y=1.05)

    # change the axis labels
    n_var=len(plot_var)
    for i in range(n_var):
        for j in range(i+1):
            #find x/y label
            xlabel = pp.axes[i][j].get_xlabel()
            ylabel = pp.axes[i][j].get_ylabel()

            # replace if found in labels dict
            if xlabel in axis_labels_dict.keys():
                pp.axes[i][j].set_xlabel(axis_labels_dict[xlabel])

            if ylabel in axis_labels_dict.keys():
                pp.axes[i][j].set_ylabel(axis_labels_dict[ylabel])

    return pp

# run plot functions
if __name__ == '__main__':

    #create parser
    parser = argparse.ArgumentParser(description='Process cuts and stuff')
    parser.add_argument("--ds", type=str, help="The dataset name (ie RD0020)",
                        required=True)
    parser.add_argument("-i", "--ion", type=str,
                        help='The ion to look at (ie "H I", "O VI")',
                        required=True)
    parser.add_argument("-r", "--refinement", type=str,
                        help="Refinement scheme, cool or nat",
                        choices=['cool', 'nat'], default='cool')
    parser.add_argument("-c", "--cut", type=str, default='cgm', nargs="*")

    args = parser.parse_args()
    #define data set info
    dsname= args.ds
    ion= args.ion
    ref = args.refinement
    cuts = args.cut

    ion_u="_".join(ion.split(" "))
    outdir=f"{homeDir}/data/absorber_data/{ref}_refinement/ion_{ion_u}/plots"

    cut_u = []
    for c in cuts:
        cut_u.append("_".join( c.split() ))

    outcut="_v_".join(cut_u)
    outdir+=f"/{outcut}"

    makedirs(outdir, exist_ok=True)
    
    #load dataframe 
    df = load_files(dsname, ion=ion, refinement=ref, cuts=cuts)
    main(df, ion, cuts, outdir)
