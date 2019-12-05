import numpy as np
import matplotlib.pyplot as plt
import yt
import trident
import sys

def main(ds_name, array_name, out_dir='./'):
    absorber_array=np.load(array_name)
    ds = yt.load(ds_fname)

    fig1 = plt.figure(figsize=(8, 6))
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

    for absorber in absorber_array:
        if absorber[3] > 15:
            curr_metal = ds.quan(absorber[7], 'dimensionless').in_units('Zsun')
            metallicity.append(np.log10(curr_metal))
            logNH.append(absorber[3])

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
    plt.scatter(slfs_logNH, slfs_metal, color='tab:grey', edgecolor='black', label='SLFSs')
    plt.scatter(plls_logNH, plls_metal, color='tab:blue', edgecolor='black', label='pLLSs')
    plt.scatter(lls_logNH, lls_metal, color='red', edgecolor='black', label='LLSs')
    plt.scatter(slls_logNH, slls_metal, color='tab:brown', edgecolor='black', label='sLLSs')
    plt.scatter(dla_logNH, dla_metal, color ='cyan', edgecolor='black', label='DLAs')
    plt.legend()
    #plt.plot(logNH, metallicity, ".")
    plt.ylabel("$log$ $\\frac{Z}{Z_{\\odot}}$", size=15)
    plt.xlabel("$log$ $N_{HI}$", size=15)
    plt.xlim(15, 22)
    plt.ylim(-3.4,0.6)

    plt.hlines([0.0, -1.4],15, 22,linestyle='dashed', alpha=0.5, zorder=1 )
    plt.title("356 Absorbers. max_impact=100kpc. z=0.3")

    #now make histogram
    fig2 = plt.figure(figsize=(8, 6))
    binwidth=0.2
    binrange=[-5, 2]
    bins = (binrange[1]-binrange[0])/binwidth
    bins = int(bins)
    plt.hist(slfs_metal, bins=bins,
             range=binrange,
             density=True,
             label=f"{len(slfs_metal)} SLFSs")
    plt.hist(plls_metal, bins=bins,
             range=binrange,
             density=True,
             linewidth=2,
             histtype='step',
             label=f"{len(plls_metal)} pLLSs")

    plt.legend()
    plt.title("Normalized Histogram of SLFSs and pLLSs")
    plt.xlabel("$log$ $\\frac{Z}{Z_{\\odot}}$", size=15)
    plt.ylabel("Normalized Count")

    fig1.save(f"{out_dir}/metallicity_H_col_density.png")
    fig2.save(f"{out_dir}/histo_SLFS_pLLS.png")

if __name__ == '__main__':
    ds_name = sys.argv[1] #"RD0036"
    max_b = sys.argv[2]
    if len(sys.argv) == 4:
        out_dir = sys.argv[3]
    else:
        out_dir='./'
    frac=0.8

    ds_fname = f"/mnt/gs18/scratch/users/boydbre1/cosmological/foggie/{dsname}/{dsname}"
    arr_dir=f"/mnt/gs18/scratch/users/boydbre1/metallicity_analysis/foggie/{dsname}/max_impact{max_b}/output{frac}"
    arr_name = f"{arr_dir}/all_absorbers.npy"

    main(ds_name, arr_name, out_dir=out_dir)
