import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import yt
import trident
import sys

def main(ds_fname, array_name, impact_param='??', out_dir='./'):
    absorber_array=np.load(array_name)
    ds = yt.load(ds_fname)

    fig = plt.figure(figsize=(12, 6))
    metallicity=[]
    logNH=[]

    slfs_logNH=[]
    plls_logNH=[]
    lls_logNH=[]
    slls_logNH=[]
    dla_logNH=[]

    slfs_metal=[]
    plls_metal=[]
    lls_metal=[]
    slls_metal=[]
    dla_metal=[]

    n_absorbers=0
    #seperate absorbers into appropriate lists
    for absorber in absorber_array:
        if absorber[3] > 15:
            curr_metal = ds.quan(absorber[7], 'dimensionless').in_units('Zsun')
            metallicity.append(np.log10(curr_metal))
            logNH.append(absorber[3])
            n_absorbers+=1

            if absorber[3] < 16.2:
                slfs_metal.append(np.log10(curr_metal))
                slfs_logNH.append(absorber[3])
            elif absorber[3] < 17.2:
                plls_metal.append(np.log10(curr_metal))
                plls_logNH.append(absorber[3])
            elif absorber[3] < 19:
                lls_metal.append(np.log10(curr_metal))
                lls_logNH.append(absorber[3])
            elif absorber[3] < 20.3:
                slls_metal.append(np.log10(curr_metal))
                slls_logNH.append(absorber[3])
            else:
                dla_metal.append(np.log10(curr_metal))
                dla_logNH.append(absorber[3])


    #plot all the pretty colors
    ax_scat = fig.add_subplot(122)
    ax_scat.scatter(slfs_logNH, slfs_metal, color='tab:grey', edgecolor='black', label='SLFSs')
    ax_scat.scatter(plls_logNH, plls_metal, color='tab:blue', edgecolor='black', label='pLLSs')
    ax_scat.scatter(lls_logNH, lls_metal, color='red', edgecolor='black', label='LLSs')
    ax_scat.scatter(slls_logNH, slls_metal, color='tab:brown', edgecolor='black', label='sLLSs')
    ax_scat.scatter(dla_logNH, dla_metal, color ='cyan', edgecolor='black', label='DLAs')
    ax_scat.legend()
    #plt.plot(logNH, metallicity, ".")
    ax_scat.set_ylabel("$log$ $\\frac{Z}{Z_{\\odot}}$", size=15)
    ax_scat.set_xlabel("$log$ $N_{HI}$", size=15)
    ax_scat.set_xlim(15, 22)
    ax_scat.set_ylim(-3.4,0.6)

    ax_scat.hlines([0.0, -1.4],15, 22,linestyle='dashed', alpha=0.5, zorder=1 )
    ax_scat.set_title(f"{n_absorbers} Absorbers. max_impact={impact_param}kpc. z={ds.current_redshift:0.2f}")

    #now make histogram
    ax_hist = fig.add_subplot(121)
    binwidth=0.2
    binrange=[-5, 2]
    bins = (binrange[1]-binrange[0])/binwidth
    bins = int(bins)
    ax_hist.hist(slfs_metal, bins=bins,
             range=binrange,
             density=True,
             label=f"{len(slfs_metal)} SLFSs")
    ax_hist.hist(plls_metal, bins=bins,
             range=binrange,
             density=True,
             linewidth=2,
             histtype='step',
             label=f"{len(plls_metal)} pLLSs")

    ax_hist.legend()
    ax_hist.set_title("Normalized Histogram of SLFSs and pLLSs")
    ax_hist.set_xlabel("$log$ $\\frac{Z}{Z_{\\odot}}$", size=15)
    ax_hist.set_ylabel("Normalized Count")

    fig.tight_layout()
    fig.savefig(f"{out_dir}/metallicity_H_col_density.png")

if __name__ == '__main__':
    ds_name = sys.argv[1] #"RD0036"
    max_b = sys.argv[2]
    if len(sys.argv) == 4:
        out_dir = sys.argv[3]
    else:
        out_dir='./'
    frac=0.8

    ds_fname = f"/mnt/gs18/scratch/users/boydbre1/cosmological/foggie/{ds_name}/{ds_name}"
    arr_dir=f"/mnt/gs18/scratch/users/boydbre1/metallicity_analysis/foggie/{ds_name}/max_impact{max_b}/output{frac}"
    arr_name = f"{arr_dir}/all_absorbers.npy"

    main(ds_fname, arr_name, impact_param=max_b, out_dir=out_dir)
