#!/bin/bash

#SBATCH --mem-per-cpu=5120mb
#SBATCH -t 12:00:00
#SBATCH -n1

s="/mnt/gs18/scratch/users/boydbre1"
dataFile="$s/isolated_galaxy/DD0076/DD0076"
mainDir="$s/multiplot_movie/movie_200kpc"
rayDir="$mainDir/rays"
ionName=$1
#label outdir but split ion into element and number
ionLabel=${ionName% *}_${ionName#* }
outDir="$mainDir/movie_${ionLabel}_images"

if ! [ -d $outDir ];
then
        mkdir -p $outDir
else
	rm -rf $outDir/*
fi

srun 									\
	-o ${outDir}/cap_movie.out 					\
	-J $ionLabel							\
	python ~/Repo/CGM/plotting_ray/capture_movie_frames.py		\
		$dataFile $rayDir "$ionName" $outDir
