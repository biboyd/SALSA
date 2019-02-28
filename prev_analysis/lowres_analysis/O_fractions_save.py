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

def _total_element_mass(field, data):
    final = data["%s_p0_mass" % element]
    index = 1
    while ("%s_p%s_mass" % (element,index)) in dir(data.ds.fields.gas):
        final += data["%s_p%s_mass" % (element,index)]
        index += 1
    return final


# By using wildcards such as ? and * with the load command, we can load up a
# Time Series containing all of these datasets simultaneously.
# '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD????/DD????'
filenames = sys.argv[1]
ts = yt.DatasetSeries(filenames,parallel=int(sys.argv[3]))

storage = {}

#use 'O VI' format
element = "O"
write_cond = sys.argv[2]

output_file="tmb_%s_fractions.dat" % (element)

if os.path.exists(output_file) and write_cond=='a':
    vals = np.loadtxt(output_file, dtype={'names': ('ds', 'time', "O_p3_mass", "O_p4_mass", "O_p5_mass", "O_p6_mass", "O_p7_mass", 'ion_mass'), 'formats': ('S6', np.float, np.float, np.float, np.float, np.float, np.float, np.float)}, 
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
    if str(ds) not in datasets and str(ds)!='DD2142':
        trident.add_ion_fields(ds, ions=[element])

        ds.add_field(('gas','total_element_mass'), 
            function=_total_element_mass, 
            units=ds.unit_system["mass"], 
            take_log=False,
            display_name='Total Element Mass')

        # Create a sphere of radius 100 kpc at the center of the dataset volume
        # sphere = ds.sphere("c", (100., "kpc"))

        ad = ds.all_data()

        element_mass = ad.quantities.total_quantity(["total_element_mass"])
        O_p3_mass = ad.quantities.total_quantity(["O_p3_mass"])
        O_p4_mass = ad.quantities.total_quantity(["O_p4_mass"])
        O_p5_mass = ad.quantities.total_quantity(["O_p5_mass"])
        O_p6_mass = ad.quantities.total_quantity(["O_p6_mass"])
        O_p7_mass = ad.quantities.total_quantity(["O_p7_mass"])


        # Calculate the entropy within that sphere
        # entr = sphere["entropy"].sum()
        # Store the current time and sphere entropy for this dataset in our
        # storage dictionary as a tuple
        store.result = (ds.current_time.in_units('Gyr'), element_mass)
        temp_time = time.time()
        if yt.is_root():
            f = open(output_file,'a')
            f.write("%s %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n" % 
                (ds, ds.current_time.in_units('Gyr'), O_p3_mass, O_p4_mass, 
                 O_p5_mass, O_p6_mass, O_p7_mass, element_mass))
            f.close()
            print("%s dataset read, finished at %f" % (ds,temp_time-start))
    
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
