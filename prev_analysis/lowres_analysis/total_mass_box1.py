import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
#plt.ioff()
import numpy as np
import yt
import trident
import time

start = time.time()
# Enable parallelism in the script (assuming it was called with
# `mpirun -np <n_procs>` )
yt.enable_parallelism()

# By using wildcards such as ? and * with the load command, we can load up a
# Time Series containing all of these datasets simultaneously.
ts = yt.load('/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD????/DD????')

storage = {}

# By using the piter() function, we can iterate on every dataset in
# the TimeSeries object.  By using the storage keyword, we can populate
# a dictionary where the dataset is the key, and sto.result is the value
# for later use when the loop is complete.

# The serial equivalent of piter() here is just "for ds in ts:" .

for store, ds in ts.piter(storage=storage):

    trident.add_ion_fields(ds, ions=['O VI'])

    # Create a sphere of radius 100 kpc at the center of the dataset volume
    # sphere = ds.sphere("c", (100., "kpc"))

    ad = ds.all_data()

    ion_mass = ad.quantities.total_quantity(["O_p5_mass"])
    #ion_mass = ad["O_p5_mass"].sum()
    # Calculate the entropy within that sphere
    # entr = sphere["entropy"].sum()
    # Store the current time and sphere entropy for this dataset in our
    # storage dictionary as a tuple
    store.result = (ds.current_time.in_units('Gyr'), ion_mass)
    temp_time = time.time()
    print("%s dataset read, at %f" % (ds,(temp_time-start)))

# Convert the storage dictionary values to a Nx2 array, so the can be easily
# plotted
if yt.is_root():
    print("did I get here yet?")
    arr = np.array(list(storage.values()))
    print("have I made this array yet?")
    print(arr)
    print(arr[:,0],arr[:,1])
    # Plot up the results: time versus entropy
    plt.semilogy(arr[:,0], arr[:,1], 'r-')
    print("have I plotted yet?")
    plt.xlabel("Time (Gyr)")
    plt.ylabel("Ion Mass (kg?)")
    plt.savefig("ion_mass_vs_time_tunedCGM.png")
    end = time.time()
    print('End at %f' % (end-start))
