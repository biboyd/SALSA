#!/bin/bash

#SBATCH --mem-per-cpu=5120mb
#SBATCH -t 12:00:00
#SBATCH -n1

scratch="/mnt/gs18/scratch/users/boydbre1"
inFile="$scratch/isolated_galaxy/DD0076/DD0076"
nRays=$1
rayDist=$2
outDir="$scratch/multiplot_movie/movie_${rayDist}kpc/rays"

if ! [ -d $outDir ]
then
	mkdir -p $outDir
else
	rm -f $outDir/*
fi

srun 								\
	-o ${outDir}/construct_rays.out 			\
	python ~/Repo/CGM/job_scripts/construct_rays.py	\
		$inFile $nRays $rayDist $outDir
