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

input_file="tmb_%s_fractions_sphere.dat" % (element)

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

vals = np.loadtxt(input_file, 
        dtype={'names': ('ds', 'time', "O_p3_mass1", "O_p4_mass1", "O_p5_mass1", "O_p6_mass1", "O_p7_mass1", 'ion_mass1',
            "O_p3_mass2", "O_p4_mass2", "O_p5_mass2", "O_p6_mass2", "O_p7_mass2", 'ion_mass2',
            "O_p3_mass3", "O_p4_mass3", "O_p5_mass3", "O_p6_mass3", "O_p7_mass3", 'ion_mass3',
            "O_p3_mass_tot", "O_p4_mass_tot", "O_p5_mass_tot", "O_p6_mass_tot", "O_p7_mass_tot", 'ion_mass_tot'), 
        'formats': ('S6', np.float, np.float, np.float, np.float, np.float, np.float, np.float,
            np.float, np.float, np.float, np.float, np.float, np.float,
            np.float, np.float, np.float, np.float, np.float, np.float,
            np.float, np.float, np.float, np.float, np.float, np.float)}, 
    comments="#", delimiter=' ')

vals=sorted(vals,key=lambda x: x[0])
print(vals)
time = [x[1]*1000 for x in vals]
O_p3_frac1 = [x[2]/x[7] for x in vals]
O_p4_frac1 = [x[3]/x[7] for x in vals]
O_p5_frac1 = [x[4]/x[7] for x in vals]
O_p6_frac1 = [x[5]/x[7] for x in vals]
O_p7_frac1 = [x[6]/x[7] for x in vals]

O_p3_frac2 = [x[8]/x[13] for x in vals]
O_p4_frac2 = [x[9]/x[13] for x in vals]
O_p5_frac2 = [x[10]/x[13] for x in vals]
O_p6_frac2 = [x[11]/x[13] for x in vals]
O_p7_frac2 = [x[12]/x[13] for x in vals]

O_p3_frac3 = [x[14]/x[19] for x in vals]
O_p4_frac3 = [x[15]/x[19] for x in vals]
O_p5_frac3 = [x[16]/x[19] for x in vals]
O_p6_frac3 = [x[17]/x[19] for x in vals]
O_p7_frac3 = [x[18]/x[19] for x in vals]

O_p3_frac_tot = [x[20]/x[25] for x in vals]
O_p4_frac_tot = [x[21]/x[25] for x in vals]
O_p5_frac_tot = [x[22]/x[25] for x in vals]
O_p6_frac_tot = [x[23]/x[25] for x in vals]
O_p7_frac_tot = [x[24]/x[25] for x in vals]
#element_mass = [x[7] for x in vals]
#print(time,element_mass)
####################
fig = plt.figure()

plt.figure(figsize=(12, 4.5))

line_3, =plt.semilogy(time, O_p3_frac1, lw=2, label='O IV')
line_4, =plt.semilogy(time, O_p4_frac1, lw=2, label='O V')
line_5, =plt.semilogy(time, O_p5_frac1, lw=2, label='O VI')
line_6, =plt.semilogy(time, O_p6_frac1, lw=2, label='O VII')
line_7, =plt.semilogy(time, O_p7_frac1, lw=2, label='O VIII')

plt.legend(handles=[line_3, line_4, line_5, line_6, line_7])

plt.ylabel("O Ion Mass Fraction [kg]")

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)

if sys.argv[2]!="none":
    y_lim = (float(sys.argv[1]),float(sys.argv[2]))
    plt.ylim(y_lim)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

plt.savefig("tmb_isogal_%s_fractions_sphere1.png" % element)
plt.close(fig)
####################
fig = plt.figure()

plt.figure(figsize=(12, 4.5))

line_3, =plt.semilogy(time, O_p3_frac2, lw=2, label='O IV')
line_4, =plt.semilogy(time, O_p4_frac2, lw=2, label='O V')
line_5, =plt.semilogy(time, O_p5_frac2, lw=2, label='O VI')
line_6, =plt.semilogy(time, O_p6_frac2, lw=2, label='O VII')
line_7, =plt.semilogy(time, O_p7_frac2, lw=2, label='O VIII')

plt.legend(handles=[line_3, line_4, line_5, line_6, line_7])

plt.ylabel("O Ion Mass Fraction [kg]")

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)

if sys.argv[2]!="none":
    y_lim = (float(sys.argv[1]),float(sys.argv[2]))
    plt.ylim(y_lim)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

plt.savefig("tmb_isogal_%s_fractions_sphere2.png" % element)
plt.close(fig)
####################
fig = plt.figure()

plt.figure(figsize=(12, 4.5))

line_3, =plt.semilogy(time, O_p3_frac3, lw=2, label='O IV')
line_4, =plt.semilogy(time, O_p4_frac3, lw=2, label='O V')
line_5, =plt.semilogy(time, O_p5_frac3, lw=2, label='O VI')
line_6, =plt.semilogy(time, O_p6_frac3, lw=2, label='O VII')
line_7, =plt.semilogy(time, O_p7_frac3, lw=2, label='O VIII')

plt.legend(handles=[line_3, line_4, line_5, line_6, line_7])

plt.ylabel("O Ion Mass Fraction [kg]")

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)

if sys.argv[2]!="none":
    y_lim = (float(sys.argv[1]),float(sys.argv[2]))
    plt.ylim(y_lim)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

plt.savefig("tmb_isogal_%s_fractions_sphere3.png" % element)
plt.close(fig)
####################
fig = plt.figure()

plt.figure(figsize=(12, 4.5))

line_3, =plt.semilogy(time, O_p3_frac_tot, lw=2, label='O IV')
line_4, =plt.semilogy(time, O_p4_frac_tot, lw=2, label='O V')
line_5, =plt.semilogy(time, O_p5_frac_tot, lw=2, label='O VI')
line_6, =plt.semilogy(time, O_p6_frac_tot, lw=2, label='O VII')
line_7, =plt.semilogy(time, O_p7_frac_tot, lw=2, label='O VIII')

plt.legend(handles=[line_3, line_4, line_5, line_6, line_7])

plt.ylabel("O Ion Mass Fraction [kg]")

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)

if sys.argv[2]!="none":
    y_lim = (float(sys.argv[1]),float(sys.argv[2]))
    plt.ylim(y_lim)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

plt.savefig("tmb_isogal_%s_fractions_sphere_tot.png" % element)
plt.close(fig)

