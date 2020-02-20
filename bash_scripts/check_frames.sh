#!/bin/bash

#runs all make movie scripts
mainDir=$1
cd $mainDir

shopt -s nullglob
for dir in movie_*_frames;
do
	num_frames=0
	for f in $dir/*.png
	do
		num_frames=$(($num_frames + 1))
	done
	echo ${dir#frames/} has  $num_frames	
done

