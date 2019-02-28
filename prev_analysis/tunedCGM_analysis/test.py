import numpy as np
import yt
import trident
yt.enable_parallelism()

print('test')

ds = yt.load('/mnt/scratch/dsilvia/simulations/reu_sims/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD1000/DD1000')
