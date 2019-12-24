# -*- coding: utf-8 -*-

"""
This is an instance of the Particle Analysis Code (PAC)

Outline: 
    1.Read tiff file sequence into aggregate objects and,
        1a. Normalize grey scale to [0,1]
        1b. Split the scan data (regions) into subregions
    2.Filter grey level intensity (GLI) data
    3.Segment GLI data using:  
        3a.Traditional topographical watershed segmentation
        3b.Power segmentation
        3c.Level-set segmentation
    4.Compute particle size and morphology distribution
    5.Locate interparticle contact
    6.Determine contact normal distribution

Classes:
    1. Particle
        Contains locations of particle voxels
        Contains morphology and size parameters

    2. Aggregate
        contains GLI data of entire region
        Contains list of all particles
        Contains grain size distribution
        Contains morphology distribution
        Contains contact normals 
        Contains list of contacting particles

    3. Reader
        Reads data from pacInput
        Crops data
        Normalizes GLI of region
        Splits region into sub-regions
        Creats aggregate objects for regions and sub-regions

    4. Filter
        Filters subregions
        Checks adequacy of filtration

    5. Segment
        Binarizes the GLI for subregions
        Checks binarizations
        Segments the particles using: 
            a. Topological watershed
            b. Power watershed
            c. Level set
        Assigns surfaces to particles

    6. Measure
        Measures particle size parameters
        Measures breakage parameters
        Measures particle morphology paramters
        Measures contact normals and Fabric
        Generates grain size distributions 
        Generates morphology distributions
        Generates contact distributions

    7. LemmeC
        Visualizes particles
        Visualizes the aggregate
        Generates plots for:
            Grain size distribution
            Contact normals
            Morphology distributions

    8. Writer
        Outputs files to pacOutput

    9. Jeeves
        Takes care of misc activites

"""

# %% Imports

print("\n\nStarting Particle analysis code (PAC)------------*")

# General imports
import numpy as np
import skimage.external.tifffile as tiffy
from skimage import io
import matplotlib.pyplot as plt
import spam.datasets as sdata

# Noun imports
import Aggregate

# Verb imports
import Reader
import Filter
import Segment
import Measure
import LemmeC
import Writer
import Jeeves

# %% Objects

reader = Reader.Reader()
filters = Filter.Filter()
segment = Segment.Segment()
measure = Measure.Measure()
lemmec = LemmeC.LemmeC()
writer = Writer.Writer()
jeeves = Jeeves.Jeeves()

# %% Read files into regions

# Folder name and image type:
folderName = '/home/eg/github/PAC/pacInput/ConeScans/OTC_90_top'
fileType = '.tiff'

# File details
pixelSize = 0.012390           # mm/pixel




# TODO: update filename to read from local folder
# TODO: Add subroutine to make grey level intensity between 0 and 1
# TODO: Read and create aggregate object, each region is a aggregate
# TODO: Subdivide region into subregions, each of which is an aggregate objects


