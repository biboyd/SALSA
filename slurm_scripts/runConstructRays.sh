#!/bin/bash

#SBATCH --mem-per-cpu=5120mb
#SBATCH -t 12:00:00
#SBATCH -n1

inFile="/mnt/home/boydbre1/data/DD0076/DD0076"
nRays=500
outDir="/mnt/home/boydbre1/data/rays"

if ! [ -d $outDir ]
then
	mkdir $outDir
else
	rm -f $outDir/*
fi

srun 								\
	-o ${outDir}/construct_rays.out 			\
	python ~/Repo/CGM/plotting_ray/test_ray_construct.py 	\
		$inFile $nRays $outDir
