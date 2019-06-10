#!/bin/bash

#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH -t 00:10:00
#SBATCH --x11=all
#SBATCH --mail-user=boydbre1@msu.edu
#SBATCH --mail-type=ALL

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

mpirun -np 16 								\
	python ~/Repo/CGM/plotting_ray/capture_movie_frames.py		\
		$dataFile $rayDir "$ionName" $outDir
scontrol show job $SLURM_JOB_ID
