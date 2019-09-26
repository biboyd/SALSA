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

#pad and annotate name

ffmpeg -i combined_ions.mp4 -vf "pad=width=4154:height=1894:x=40:y=75:color=white" padded_ions.mp4

#ffmpeg -i padded_ions.mp4 -vf drawtext="fontsize=30:font=/usr/share/fonts/truetype/freefont/FreeMono.ttf:text='FOGGIE - B. Boyd(MSU)':x=(w-text_w)/1.01:y=(h-text_h)/1.01" final_ions.mp4

