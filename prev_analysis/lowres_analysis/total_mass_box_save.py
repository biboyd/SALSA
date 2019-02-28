import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

import yt
import trident
import time
import os, sys

# Enable parallelism in the script (assuming it was called with
# `mpirun -np <n_procs>` )
start = time.time()

yt.enable_parallelism()

# return formatted string for ion input string
def ion_p_name(ion):
    return ion.split()[0]+'_p'+str(trident.from_roman(ion.split()[1])-1)

# remove space from yt quan; accepts h5file
def yt_remove_spaces(ytquan):
    return str(ytquan).replace(" ", "")

# By using wildcards such as ? and * with the load command, we can load up a
# Time Series containing all of these datasets simultaneously.
# '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD????/DD????'
filenames = sys.argv[1]
ts = yt.DatasetSeries(filenames,parallel=int(sys.argv[4]))

storage = {}
#use 'O VI' format
ion = sys.argv[2]
write_cond = sys.argv[3]
#filename = open("ERRORFILES_%s.txt" % ion_p_name(ion),'w')
#filename.close()
if ion == "none":
    output_file="tmb.dat"
    ion = "O VI"
else:
    output_file="tmb_%s.dat" % ion_p_name(ion)

if os.path.exists(output_file) and write_cond=='a':
    vals = np.loadtxt(output_file, dtype={'names': ('ds', 'time', 'ion_mass'), 'formats': ('S6', np.float, np.float)}, 
    comments="#", delimiter=' ')
    vals=sorted(vals,key=lambda x: x[0])
    datasets = [x[0].decode('utf-8') for x in vals]
elif write_cond=='w':
    if yt.is_root():
        f = open(output_file,'w')
        f.close()
    datasets = []
else:
    print('Error -- use correct write condition')
    sys.exit()


for store, ds in ts.piter(storage=storage):
    if str(ds) not in datasets:
        print("Working on",ds)
        trident.add_ion_fields(ds, ions=[ion])

        # Create a sphere of radius 100 kpc at the center of the dataset volume
        # sphere = ds.sphere("c", (100., "kpc"))

        ad = ds.all_data()

        #try:
        if True:
            ion_mass = ad.quantities.total_quantity(["%s_mass" % ion_p_name(ion)])

        # Calculate the entropy within that sphere
        # entr = sphere["entropy"].sum()
        # Store the current time and sphere entropy for this dataset in our
        # storage dictionary as a tuple
            store.result = (ds.current_time.in_units('Gyr'), ion_mass)
            temp_time = time.time()
            if yt.is_root():
                f = open(output_file,'a')
                f.write("%s %.6e %.6e\n" % (ds, ds.current_time.in_units('Gyr'), ion_mass))
                f.close()
                print("%s dataset read, finished at %f" % (ds,temp_time-start))
        #except:
        #    print("%s FAILED" % ds)
        #    if yt.is_root():
        #        filename = open("ERRORFILES_%s.txt" % ion_p_name(ion),'w')
        #        filename.write("%s failed\n" % ds)  
        #        filename.close()
        #    print("SHOULD HAVE WRITTEN\n\n\n")
        #    sys.exit()
        #print("test 1")
    #print("test 2")
#print("test 3")
#quit()
#sys.exit()
#print("test 4")
# Convert the storage dictionary values to a Nx2 array, so the can be easily
# plotted
#if yt.is_root():
#    arr = np.array(list(storage.values()))
#    print(arr)
#    # Plot up the results: time versus entropy
#    plt.semilogy(arr[:,0], arr[:,1], 'r-')
#    plt.xlabel("Time (Gyr)")
#    plt.ylabel("Ion Mass (kg?)")
#    plt.savefig("ion_mass_vs_time_O_p5_lowres.png")
#    end = time.time()
#    print("Completely finished at %f s" %(end-start))
