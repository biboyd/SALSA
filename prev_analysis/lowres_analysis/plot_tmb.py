import matplotlib as mpl
mpl.use("agg")

import yt
import trident
import numpy as np
import sys
import matplotlib.pyplot as plt

# return formatted string for ion input string
def ion_p_name(ion):
    return ion.split()[0]+'_p'+str(trident.from_roman(ion.split()[1])-1)

ion = sys.argv[1]

if ion == "none":
    input_file="tmb.dat"
else:
    input_file="tmb_%s.dat" % ion_p_name(ion)

mpl.rcParams['font.family'] = 'serif'
mpl.rc('xtick', labelsize=15)
mpl.rc('ytick', labelsize=15)
mpl.rc('axes', labelsize=18)
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['legend.fontsize'] = 18
mpl.rcParams['xtick.major.width']= 1.0
mpl.rcParams['xtick.minor.width']= 1.0
mpl.rcParams['xtick.major.size']= 8.0
mpl.rcParams['xtick.minor.size']= 4.0
mpl.rcParams['ytick.major.width']= 1.0
mpl.rcParams['ytick.minor.width']= 1.0
mpl.rcParams['ytick.major.size']= 8.0
mpl.rcParams['ytick.minor.size']= 4.0

vals = np.loadtxt(input_file, dtype={'names': ('ds', 'time', 'ion_mass'), 'formats': ('S6', np.float, np.float)}, 
    comments="#", delimiter=' ')

vals=sorted(vals,key=lambda x: x[0])
print(vals)
time = [x[1]*1000 for x in vals]
ion_mass = [x[2] for x in vals]
print(time,ion_mass)

plt.figure(figsize=(12, 4.5))

plt.plot(time, ion_mass, lw=2,color='k')
#plt.semilogy(time, ion_mass, lw=2, color='k')

if ion == "none":
    plt.ylabel("Total O VI Mass [kg]")
else:
    plt.ylabel("Total %s Mass [kg]" % ion)

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)

if sys.argv[2]!="none":
    y_lim = (float(sys.argv[2]),float(sys.argv[3]))
    plt.ylim(y_lim)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

if ion == "none":
	plt.savefig("tmb_isogal.png")
else:
    plt.savefig("tmb_isogal_%s.png" % ion_p_name(ion))
