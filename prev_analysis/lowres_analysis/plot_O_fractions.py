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

element = 'O'

input_file="tmb_%s_fractions.dat" % (element)

mpl.rcParams['font.family'] = 'serif'
mpl.rc('xtick', labelsize=15)
mpl.rc('ytick', labelsize=15)
mpl.rc('axes', labelsize=18)
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['xtick.major.width']= 1.0
mpl.rcParams['xtick.minor.width']= 1.0
mpl.rcParams['xtick.major.size']= 8.0
mpl.rcParams['xtick.minor.size']= 4.0
mpl.rcParams['ytick.major.width']= 1.0
mpl.rcParams['ytick.minor.width']= 1.0
mpl.rcParams['ytick.major.size']= 8.0
mpl.rcParams['ytick.minor.size']= 4.0

vals = np.loadtxt(input_file, dtype={'names': ('ds', 'time', "O_p3_mass", "O_p4_mass", "O_p5_mass", "O_p6_mass", "O_p7_mass", 'ion_mass'), 'formats': ('S6', np.float, np.float, np.float, np.float, np.float, np.float, np.float)}, 
    comments="#", delimiter=' ')

vals=sorted(vals,key=lambda x: x[0])
print(vals)
time = [x[1]*1000 for x in vals]
O_p3_frac = [x[2]/x[7]*10 for x in vals]
O_p4_frac = [x[3]/x[7]*10 for x in vals]
O_p5_frac = [x[4]/x[7]/1 for x in vals]
O_p6_frac = [x[5]/x[7]/100 for x in vals]
O_p7_frac = [x[6]/x[7]/100 for x in vals]
#element_mass = [x[7] for x in vals]
#print(time,element_mass)

plt.figure(figsize=(12, 4.5))

#line_3, = plt.semilogy(time, O_p3_frac, lw=2, label='O IV')
#line_4, =plt.semilogy(time, O_p4_frac, lw=2, label='O V')
#line_5, =plt.semilogy(time, O_p5_frac, lw=2, label='O VI')#line_6, =plt.semilogy(time, O_p6_frac, lw=2, label='O VII')
#line_6, =plt.semilogy(time, O_p6_frac, lw=2, label='O VII')
#line_7, =plt.semilogy(time, O_p7_frac, lw=2, label='O VIII')

line_3, = plt.plot(time, O_p3_frac, lw=2, label='O IV * 10')
line_4, =plt.plot(time, O_p4_frac, lw=2, label='O V * 10')
line_5, =plt.plot(time, O_p5_frac, lw=2, label='O VI')
line_6, =plt.plot(time, O_p6_frac, lw=2, label='O VII / 100')
line_7, =plt.plot(time, O_p7_frac, lw=2, label='O VIII / 100')


plt.legend(loc=2, handles=[line_3, line_4, line_5, line_6, line_7])

plt.ylabel("O Ion Mass Fraction [kg]")

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)

if sys.argv[2]!="none":
    y_lim = (float(sys.argv[1]),float(sys.argv[2]))
    plt.ylim(y_lim)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

plt.savefig("tmb_isogal_%s_fractions_scaled.png" % element)
