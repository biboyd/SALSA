#!/bin/bash

#runs all make movie scripts
frameDir=$1
force="$2"
cd $frameDir

#check for equal number of frames
H_frames=`ls movie_H_I_frames/*.png -1 | wc -l`
C_frames=`ls movie_C_IV_frames/*.png -1 | wc -l`
O_frames=`ls movie_O_VI_frames/*.png -1 | wc -l`

if [ $H_frames == $C_frames ] && [ $C_frames == $O_frames ]
then
	echo "total of $H_frames each"
else
	echo "Number of frames must be equal"
	echo "H I  has $H_frames"
	echo "C IV has $C_frames"
	echo "O VI has $O_frames"
        #check if force movie was called (-f)
        if [ "$force" != "-f" ]
	then
		exit 1
	fi
fi
#make individual movies using movie.sh
for i in H_I C_IV O_VI Si_II Si_III Ne_VIII Mg_X
do 
	~/Repo/CGM/bash_scripts/movie.sh movie_${i}.mp4 5 movie_${i}_frames/*.png 
done

#for i in H_I C_IV O_VI
#do 
#        ~/Repo/CGM/movie_scripts/movie.sh movie_${i}_slow.mp4 1 frames/movie_${i}_frames/*.png 
#done

#make combined movies

ffmpeg  -y						 \
	-i movie_H_I.mp4 				 \
	-i movie_C_IV.mp4				 \
	-i movie_O_VI.mp4				 \
	-filter_complex "[0:v][1:v][2:v]hstack=inputs=3" \
	./combined_HCO_movie.mp4

#ffmpeg  -y 						 \
#	-i movie_H_I_slow.mp4 				 \
#	-i movie_C_IV_slow.mp4 				 \
#	-i movie_O_VI_slow.mp4 				 \
#	-filter_complex "[0:v][1:v][2:v]hstack=inputs=3" \
#	./combined_movie_slow.mp4

