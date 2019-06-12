# Code Related to analyzing the CGM by creating Synthetic Spectra
## using enzo, yt, and trident projects

# Directory layout
## prev_analysis
Holds code carried over from a previous analysis
some code is borrowed/modified

## impact_par_selection 
Holds code related to verifying that the impact parameter for sightlines are being sampled in an appropriate manner that mirrors observations.

## plotting_ray
Holds code for a multiplot class that uses yt and trident. It's purpose is to construct a multiplot that shows a light ray traveling through a galactic simulation, it's corresponding number density along its length (for a given ion), and the resulting spectra an observer would see.  
main file is plotter.py (holds relavent classes). 
construct_rays.py holds code for making the rays used in the movie making
and capture_movie_frames.py hold code for running the jobs to make the movie. 


## misty_fog_spectacle
Holds (currently deprecated) wrapper for extracting information from trident data including converting from 
wavelength space to velocity space. 
Additionally there is code for experimenting with Spectacle (spectra fitting software)

## slurm_scripts
scripts for submitting jobs related to making a movie out of the multiplot
