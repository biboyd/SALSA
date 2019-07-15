#!/bin/bash

mainDir=$1

cd $mainDir


#stack cold/cool warm/hot horizontal
ffmpeg  -y						\
	-i movie_cold_gas.mp4 				\
	-i movie_cool_gas.mp4				\
	-filter_complex "[0:v][1:v]hstack=inputs=2"     \
	./tmp_c.mp4

ffmpeg  -y						\
	-i movie_warm_gas.mp4 				\
	-i movie_hot_gas.mp4 				\
	-filter_complex "[0:v][1:v]hstack=inputs=2"     \
	./tmp_h.mp4

#stack vertically
ffmpeg  -y						\
	-i tmp_c.mp4                                    \
	-i tmp_h.mp4					\
	-filter_complex "[0:v][1:v]vstack=inputs=2"     \
	./combined_gas.mp4

rm tmp_?.mp4
