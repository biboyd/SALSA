#!/bin/bash

#SBATCH --mem-per-cpu=5120mb
#SBATCH -t 12:00:00
#SBATCH -n1

inFile="/mnt/home/boydbre1/data/DD0076/DD0076"
nRays=$1
outDir="/mnt/home/boydbre1/data/Closer_rays"

if ! [ -d $outDir ]
then
	mkdir $outDir
else
	rm -f $outDir/*
fi

srun 								\
	-o ${outDir}/construct_rays.out 			\
	python ~/Repo/CGM/plotting_ray/construct_rays.py	\
		$inFile $nRays $outDir
