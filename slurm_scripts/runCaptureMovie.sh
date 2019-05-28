#!/bin/bash

#SBATCH --mem-per-cpu=5120mb
#SBATCH -t 12:00:00
#SBATCH -n1

tilda="/mnt/home/boydbre1"
dataFile="$tilda/data/DD0076/DD0076"
inDir="$tilda/data/rays"
ionName=$1
#label outdir but split ion into element and number
ionLabel=${ionName% *}_${ionName#* }
outDir="$tilda/data/movie_${ionLabel}_images"

if ! [ -d $outDir ];
then
        mkdir $outDir
fi

srun 									\
	-o ${outDir}/cap_movie.out 					\
	--x11=all							\
	python ~/Repo/CGM/plotting_ray/test_movie_class.py 		\
		$dataFile $inDir "$ionName" $outDir
