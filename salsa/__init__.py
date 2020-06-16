"""
salsa
Written by: Brendan Boyd boyd.brendan@stonybrook.edu

"""

__version__ = "0.0.0"

from salsa.absorber_extractor import absorber_extractor

from salsa.absorber_plotter import absorber_plotter

from salsa.generate_light_rays import generate_lrays, construct_rays, random_sightlines

from salsa.generate_catalog import generate_catalog

from salsa.utils.filter_definitions import parse_cut_filter
