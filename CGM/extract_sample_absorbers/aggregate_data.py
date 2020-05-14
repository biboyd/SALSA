from sys import argv
from os import makedirs
from CGM.general_utils.collect_files import combine_astropy_files

inDir = argv[1]
ds = argv[2]
cuts = argv[3]
outDir = argv[4]

cut_dir = "_".join( cuts.split(" ") )

inDir+= f"/{cut_dir}"
outDir += f"/{cut_dir}"
makedirs(outDir, exist_ok=True)

outfile=f"{outDir}/{ds}_absorbers.h5"
combine_astropy_files(inDir, outfile=outfile)
