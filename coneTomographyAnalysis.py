'''
Code for the analysis of 3D tomo data for samples collected around the cone penetrometer

written by @eg
'''

import numpy as np
import skimage.external.tifffile as tf
import time

from pac import Reader
from pac import Segment
from pac import Measure
from pac import Plot

# Each test has three samples - top mid and bottom
scanData = ['2QR25', '2QR50', '2QR90','OGF25', 'OGF50', 'OGF90','OTC25', 'OTC50', 'OTC90']

for sample in scanData:
    '''
    Initialize the samples variables
    Read the files
    Extract subregions from the 3D tomo data - 6D50 edge length
    Compute particle sizes of the sand samples - user threshold after initial guess from OTSU
    Get particel aspect ration from PCA length ratios
    Get particle size distribution
    Calculate relative breakage parameter
    ---
    Get contact normals
    Get fabric tensors
    Get fabric projection plots
    '''

    # Initialize sample variables
    #   Input file location (3 for top, mid, and tip)
    #   Output file location
    #   Resolution
    #   Subregion size
    #   Shouler reference location in the top scan - this is a point which allows measurement of distances
    #   Starting subregion locations in each scan

    # Read the files
    #   Read one scan at a time and only the center region large enough to have the subregion size
    #   Extract subregions in a loop till all the subregions needed are analyzed. - no fabric so the gradation analysis should be faster

    # Get particle size distribution for subregions
    #   Take subregions from the file and analyze for particle size distribution
    #   No physical threhold - obtain intial guess form OTSU for threshold and then iterate till:
    #       1. Capture the largest particle size (using Feret diameter) same as the original gradation
    #       2. Visally no errors in segmentation - over or under segmentation

    # Get from particle size distribution:
    #   Relative breakage parameter
    #   Morphology distribution

    # Get from labelled map:
    #   Contact normal distribution
    #   Fabric tensors
    #   Coordination number distribution




