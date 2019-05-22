# Code Related to analyzing the CGM
## using enzo, yt, and trident projects

# Directory layout
## prev_analysis
Holds code carried over from a previous analysis
some code is borrowed/modified

## impact_par_selection 
Holds code related to verifying that the impact parameter for sightlines are being sampled in an appropriate manner that mirrors observations.

## plotting_ray
Holds code for a multiplot class that uses yt and trident. It's purpose is to construct a multiplot that shows a light ray traveling through a galactic simulation, it's corresponding number density along its length (for a given ion), and the resulting spectra an observer would see.
