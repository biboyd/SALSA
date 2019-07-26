#!/bin/bash

mainDir=/mnt/gs18/scratch/users/boydbre1/multiplot_movie/foggie/

shopt -s nullglob

for ds in $mainDir/RD*;
do
	echo ${ds#$mainDir/}:
	for dist in $ds/movie_*kpc;
	do
		cnt=0
		for r in $dist/rays/*.h5;
		do
			cnt=$((cnt+1))
		done
		echo ${dist#$ds/movie_} rays=$cnt 
	done
        echo ''
done
 
