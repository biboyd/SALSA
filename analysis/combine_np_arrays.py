import numpy as np
from os import listdir
from sys import argv

def combine_arrays(dir_name):
    all_files = listdir(dir_name)
    
    #files to exclude
    black_list=["absorber_info_header.npy", "all_absorbers.npy"]

    # get np files only
    np_files=[]
    for f in all_files:
        if f[-3:]=='npy' and f not in black_list:
            np_files.append(np.load(f"{dir_name}/{f}"))
    
    all_abs = np.vstack(np_files)
    
    return all_abs

if __name__ == '__main__':
    dir_name=argv[1]

    all_abs = combine_arrays(dir_name)
    np.save(f"{dir_name}/all_absorbers.npy", all_abs)

