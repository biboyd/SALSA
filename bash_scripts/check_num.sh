#!/bin/bash

#runs all make movie scripts
mainDir=$1
ext=$2
shopt -s nullglob
num_frames=0
for f in $mainDir/*$ext
do
	num_frames=$(($num_frames + 1))
done
echo $mainDir has  $num_frames	

