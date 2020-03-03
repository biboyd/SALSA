import numpy as np
from os import listdir
from sys import argv
dir_name=argv[1]

all_files = listdir(dir_name)

# get np files only
np_files=[]
for f in all_files:
    if f[-3:]=='npy' and f != "absorber_info_header.npy":
        np_files.append(np.load(f"{dir_name}/{f}"))

all_abs = np.vstack(np_files)

np.save(f"{dir_name}/all_absorbers.npy", all_abs)
