import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt
import sys

from os import makedirs
import seaborn as sns

## MAKE COLUMN DENSITY HISTO comparing Spectacle and Cloudy
def col_dense_histo(spect_df, ice_df, z,ion="O VI", bins=20, hist_range=None, outdir='./'):
    ion_u=f"{ion.split()[0]}_{ion.split()[1]}"
    spect_cd = spect_df['col density']
    ice_cd = ice_df['col density']

    fig, ax = plt.subplots(1)
    double_hist(spect_cd, ice_cd, bins, hist_range, ax=ax, color1='tab:purple', color2='black')

    ax.set_title(f"{ion} Column Density for Spectacle and ICE (z={z:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel(f"Log( $N_{ {ion} }$ )")
    ax.legend()

    fig.savefig(f"{outdir}/col_density_{ion_u}_{z:0.2f}.png", dpi=400)

## MAKE VELOCITY HISTO
def los_velocity_histo(spect_df, ice_df, z,ion="O VI", bins=25, hist_range=None,outdir='./'):
    ion_u=f"{ion.split()[0]}_{ion.split()[1]}"
    spect_vel = spect_df['avg_velocity']
    ice_vel = ice_df['avg_velocity']

    fig, ax = plt.subplots(1)
    double_hist(spect_vel, ice_vel, bins, hist_range, ax=ax, color1='tab:purple', color2='black')

    ax.set_title(f"{ion} Velocity for Spectacle and ICE (z={z:0.2f})")
    ax.set_ylabel("% of absorbers")
    ax.set_xlabel("Line of Sight Velocity (km/s)")
    ax.legend()

    fig.savefig(f"{outdir}/velocity_{ion_u}_{z:0.2f}.png", dpi=400)
