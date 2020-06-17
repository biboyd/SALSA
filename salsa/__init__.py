"""
salsa
Written by: Brendan Boyd boyd.brendan@stonybrook.edu

"""

__version__ = "0.0.0"

from salsa.absorber_extractor import AbsorberExtractor

from salsa.absorber_plotter import AbsorberPlotter

from salsa.generate_light_rays import generate_lrays, construct_rays, random_sightlines

from salsa.generate_catalog import generate_catalog, get_catalog

from salsa import utils

from salsa.utils.filter_definitions import parse_cut_filter
