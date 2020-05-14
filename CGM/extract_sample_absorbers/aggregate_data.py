from sys import argv
from CGM.general_utils.collect_files import combine_astropy_files

inDir = argv[1]
ds = argv[2]
cuts = argv[3]
outDir = argv[4]

cut_dir = "_".join( cuts.split(" ") )

inDir+= f"/{cut_dir}"
outfile = f"{outDir}/{cut_dir}/{ds}_absorbers.h5"

combine_astropy_files(inDir, outfile=outfile)
