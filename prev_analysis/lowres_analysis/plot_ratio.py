import matplotlib as mpl
mpl.use("agg")

import yt
import trident
import numpy as np
import sys
import matplotlib.pyplot as plt

# return formatted string for  input string
def _p_name(ion):
    return ion.split()[0]+'_p'+str(trident.from_roman(ion.split()[1])-1)

input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]
y_label = sys.argv[4]

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

vals1 = np.loadtxt(input_file1, dtype={'names': ('ds', 'time', 'mass1'), 'formats': ('S6', np.float, np.float)}, 
    comments="#", delimiter=' ')
vals1=sorted(vals1,key=lambda x: x[0])
#TRY THIS METHOD
seen = set()
uniq1 = []
for i in np.arange(0,len(vals1)):
    x=(vals1[i][1],vals1[i][2])
    if x not in seen:
        uniq1.append(x)
        seen.add(x)
uniq1=sorted(uniq1,key=lambda x: x[0])

time1 = [x[0]*1000 for x in uniq1]
mass1 = [x[1] for x in uniq1]
print("Sizes (1):",len(time1),len(mass1))

# for reading sfr
# vals = np.loadtxt("sfrd.dat", comments="#", delimiter=' ')

# time = vals[:,1]
# sfr = vals[:,2]

vals2 = np.loadtxt(input_file2, dtype={'names': ('ds', 'time', 'mass2'), 'formats': ('S6', np.float, np.float)}, 
    comments="#", delimiter=' ')
vals2=sorted(vals2,key=lambda x: x[0])

# seen = set()
# uniq2 = []
# for i in np.arange(0,len(vals2)):
#     x=(vals2[i][1],vals2[i][2])
#     if x not in seen:
#         uniq2.append(x)
#         seen.add(x)
# uniq2=sorted(uniq2,key=lambda x: x[0])

time2 = [x[1]*1000 for x in vals2]
mass2 = [x[2] for x in vals2]
print("Sizes (2):",len(time2),len(mass2))

# file = open("WHYWHYWHY1.txt",'a')
# for (x,y) in zip(time1,mass1):
#     file.write(str(x)+' '+str(y)+"\n")
# file.close()

# file = open("WHYWHYWHY2.txt",'a')
# for (x,y) in zip(time2,mass2):
#     file.write(str(x)+' '+str(y)+"\n")
# file.close()


if time1 == time2:
    ratio = [x/y for x, y in zip(mass1, mass2)]
else:
    print("Error -- time mismatch.")
    sys.exit()

plt.figure(figsize=(12, 4.5))

plt.semilogy(time1, ratio, lw=2, color='k')

plt.ylabel(str(y_label))

plt.xlabel("Time [Myr]")

plt.xlim(0.0, 12000.0)

if sys.argv[5]!="none":
    y_lim = (float(sys.argv[5]),float(sys.argv[6]))
    plt.ylim(y_lim)

ax = plt.axes()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.025, box.width * 1.05, box.height * 1.05])

plt.savefig(output_file)
