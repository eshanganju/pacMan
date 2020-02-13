# -*- coding: utf-8 -*-

"""
This is an instance of the Particle Analysis Code (PAC)

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


