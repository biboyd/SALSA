#!/bin/bash

mainDir=$1

cd $mainDir

#stack cold/cool warm/hot horizontal
ffmpeg  -y                                              \
        -i movie_density.mp4                           \
        -i movie_Si_p1_number_density.mp4                           \
        -i movie_Si_p2_number_density.mp4                           \
        -i movie_Si_p3_number_density.mp4                           \
        -filter_complex "[0:v][1:v][2:v][3:v]hstack=inputs=4"     \
        ./tmp_top.mp4

ffmpeg  -y                                              \
        -i movie_H_p0_number_density.mp4                           \
        -i movie_C_p1_number_density.mp4                           \
        -i movie_C_p3_number_density.mp4                           \
        -i movie_O_p5_number_density.mp4                           \
        -filter_complex "[0:v][1:v]hstack=inputs=4"     \
        ./tmp_bottom.mp4

ffmpeg	-y 						\
	-i tmp_top.mp4					\
	-i tmp_bottom.mp4				\
	-filter_complex "[0:v][1:v]vstack=inputs=2"	\
	./combined_ions.mp4

rm tmp_top.mp4 tmp_bottom.mp4
