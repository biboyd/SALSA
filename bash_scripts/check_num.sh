#!/bin/bash

#runs all make movie scripts
mainDir=$1

shopt -s nullglob
num_frames=0
for f in $mainDir/*.png
do
	num_frames=$(($num_frames + 1))
done
echo $mainDir has  $num_frames	

