"""
Synthetic Absorption Line Surveyor Application (SALSA) for extracting absorbers from
sightlines in astrophysical hydrodynamics simulations.
"""

from salsa.absorber_extractor import AbsorberExtractor, SPICEAbsorberExtractor
from salsa.absorber_plotter import AbsorberPlotter

from salsa.generate_light_rays import generate_lrays, construct_rays, random_sightlines

from salsa.generate_catalog import generate_catalog

from salsa import utils
