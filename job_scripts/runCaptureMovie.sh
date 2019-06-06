#!/bin/bash

#SBATCH --mem-per-cpu=5120mb
#SBATCH -t 12:00:00
#SBATCH -n1
#SBATCH --x11=all

ionName=$1
rayDist=$2
s="/mnt/gs18/scratch/users/boydbre1"
dataFile="$s/isolated_galaxy/DD0076/DD0076"
mainDir="$s/multiplot_movie/movie_${rayDist}kpc"
rayDir="$mainDir/rays"
#label outdir but split ion into element and number
ionLabel=${ionName% *}_${ionName#* }
outDir="$mainDir/frames/movie_${ionLabel}_frames"

if ! [ -d $outDir ];
then
        mkdir -p $outDir
else
	rm -rf $outDir/*
fi

srun 									\
	-o ${outDir}/logs.out 						\
	python ~/Repo/CGM/job_scripts/capture_movie_frames.py		\
		$dataFile $rayDir "$ionName" $outDir
#combine images to make movie
~/Repo/CGM/movie.sh movie_${ionLabel}.mp4 5 $outDir/*.png

