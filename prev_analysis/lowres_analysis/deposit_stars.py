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

@yt.particle_filter(requires=["particle_type"], filtered_type='all')
def stars(pfilter, data):
    filter = data[(pfilter.filtered_type, "particle_type")] == 2
    return filter

filenames = sys.argv[1]
ts = yt.DatasetSeries(filenames,parallel=int(sys.argv[2]))

for ds in ts.piter():
    # example load
    #ds = yt.load('/mnt/scratch/dsilvia/simulations/reu_sims/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD1000/DD1000')
    if str(ds)=='DD0002' or str(ds)=='DD0003':
        continue
    stars_exist = ds.add_particle_filter('stars')

    proj_filename = "DEPOSIT_STARS/deposit_stars_extended_z_%s.png" % ds
    if not os.path.exists(proj_filename) and stars_exist:
        p = yt.ProjectionPlot(ds,'z',('deposit','stars_cic'),width=(600,'kpc'),center='c')
        
        ad = ds.all_data()
        pix_num=1000
        res = [pix_num, pix_num]
        width = ds.quan(600,'kpc')
        frb = p.data_source.to_frb(width, res)
        proj_array = np.array(frb[('deposit','stars_cic')])
        max_density = np.max(proj_array)
        
        #max_density = max(p[('deposit','stars_cic')])
        p.annotate_title("Max Density is %.3e g/cm^2" % max_density)
        p.set_cmap(('deposit','stars_cic'),'STD GAMMA-II')
        p.set_zlim(('deposit','stars_cic'),1e-16,1e0)
        p.annotate_sphere([0.5, 0.5, 0.5], radius=(100, 'kpc'),
                          circle_args={'color':'black'})
        p.annotate_sphere([0.5, 0.5, 0.5], radius=(200, 'kpc'),
                          circle_args={'color':'black'})
        p.annotate_sphere([0.5, 0.5, 0.5], radius=(300, 'kpc'),
                          circle_args={'color':'black'})
        p.annotate_timestamp(redshift=False,draw_inset_box=True)
        p.save(proj_filename)
