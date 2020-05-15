# Code Related to analyzing the CGM by creating Synthetic Spectra
using yt, trident and spectacle projects

# Directory layout
## CGM
holds code to construct a python package that includes `absorber_extractor` and `absorber_plotter` classes that can extract absorbers from a simulation given trident lightrays. Also contains scripts for extracting distribution of absorbers and analyzing them as well as scripts seperate from that main class for visualizing simulation.  

### CGM/absorber_analysis
Scripts for creating plots of different absorption properties. Particularly for exploring as well as comparing different type of absorbers ie hot/cold or inflow/outflow.

### CGM/absorber_extraction_class
This contains the code for the `absorber_extractor` which can extract absorbers from a trident lightray using either spectacle to fit spectra or the "ICE Method" which uses an algorithm to find the absorbers directly from the number density along the lightray. 
### CGM/comparison_code 
A script comparing the results of different algorithms for extracting absorbers directly from lightray

### CGM/extract_sample_absorbers
Scripts to construct lightrays and then extract absorbers from those lightrays. One creates a uniform "scan" of the galaxy which can then be constructed into a movie. The other creates a random distribution of lightrays which then can be used to find a distribution of absorbers that can be compared to observational studies. Cuts can be made to the dataset before extracting absorbers which allows to study different type of absorbers.

### CGM/general_utils
Contains some scripts/functions that different scripts all use. Additionally `filter_definitions.py` contains different dictionaries with default values that are used by different scripts

### CGM/notebooks
notebook to explain how to use the `absorber_extractor`/`absorber_plotter` classes. (Currently very out of date)

### CGM/visualize_gal
Scripts for visualizing the galaxy. `create_slices` creates slices going through the galaxy as well as a density projection face on. `create_rot_frames` creates a series of projection of the galaxy at different angles so they can be combined into a movie that effectively rotates the galaxy and give a good 3d representation.

## multi_plot_movie
Holds code for a multiplot class that uses yt and trident. It's purpose is to construct a multiplot that shows a light ray traveling through a galactic simulation, it's corresponding number density along its length (for a given ion), and the resulting spectra an observer would see.  
* main file is multi_plot.py (holds relavent classe).   
* construct_rays.py holds code for making the rays used in the movie making  
* capture_movie_frames.py hold code for running the jobs to make the movie.   
Also holds center_finder.py which reads/records the center, orientation, rshift, and bulk velocity for each galaxy  

## job_scripts
scripts for submitting jobs to slurm   

## polar_dependence
Code for looking at the relation between column density and polar angle
from axis of galaxy. (likely deprecated) 

## bash_scripts
* scripts for making movies (combining frames) and for editing movies (combining side by side etc.)  
* scripts for checking num frames/rays already made. helpful to figure out if/what jobs failed 
