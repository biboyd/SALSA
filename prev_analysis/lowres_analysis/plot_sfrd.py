import matplotlib as mpl
mpl.use("agg")

import yt
import numpy as np

import matplotlib.pyplot as plt

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

vals = np.loadtxt("sfrd.dat", comments="#", delimiter=' ')

time = vals[:,1]
sfr = vals[:,2]

plt.figure(figsize=(12, 4.5))

plt.semilogy(time, sfr, lw=2, color='k')

plt.ylabel("SFR [Msun/yr]")

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)
plt.ylim(1e-2, 1e2)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

plt.savefig("sfr_isogal.png")
