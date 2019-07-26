# Code Related to analyzing the CGM by creating Synthetic Spectra
## using enzo, yt, trident and spectacle projects

# Directory layout
## prev_analysis
Holds code carried over from a previous analysis
some code is borrowed/modified  

## multi_plot
Holds code for a multiplot class that uses yt and trident. It's purpose is to construct a multiplot that shows a light ray traveling through a galactic simulation, it's corresponding number density along its length (for a given ion), and the resulting spectra an observer would see.  
main file is multi_plot.py (holds relavent classe).   
construct_rays.py holds code for making the rays used in the movie making  
and capture_movie_frames.py hold code for running the jobs to make the movie.   
Also holds center_finder.py which reads/records the center, orientation, rshift, and bulk velocity for each galaxy  


## spectacle_testing 
Holds (currently deprecated) wrapper for extracting information from trident data including converting from 
wavelength space to velocity space.   
Additionally there is code for experimenting with Spectacle (spectra fitting software)

## job_scripts
scripts for submitting jobs to slurm   

## phi_dependence
Code for looking at the relation between column density and polar angle
from axis of galaxy.   

## rotation_movie
Code for making movies of rotating projections of a galaxy. Useful for seeing the structure 

## bash_scripts
scripts for making movies (combining frames) and for editing movies (combining side by side etc.)  
scripts for checking num frames/rays already made  
