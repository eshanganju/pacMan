# -*- coding: utf-8 -*-

"""
This is the second version of the Particle analysis code

The following are attempted with trepidation: 
  1.Read tiff files containing CT scan data in to np arrays and
    1a.Normalize grey scale to [0,1]
  2.Filter CT scan data to reduce noise using automated NLM loops
  3.Segment CT data using:  
    3a.Traditional topographical watershed segmentation
    3b.Power segmentation
    3c.Level-set segmentation
  4.Compute particle size and morphology distribution
  5.Locate interparticle contact
  6.Determine contact normal distribution

Classes:
    1. Particle
    2. Aggregate
---
    3. Reader
    4. Filter
    5. Segment
    6. Measure
    7. LemmeC
    8. Writer
"""

# %% General Classes
import numpy as np
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt
import spam.datasets as sdata

# %% Unique Classes
import Reader
import Aggregate
import Filter
import Segment
import Measure
import LemmeC
import Writer

# %% Objects

print("\n\nStarting Particle analysis code (PAC)------------*")
reader = Reader.Reader()
filters = Filter.Filter()
segment = Segment.Segment()
measure = Measure.Measure()
lemmec = LemmeC.LemmeC()
writer = Writer.Writer()

# %% File
fileName = 'Box4B.tiff'
# update filename to read from local folder

pixelSize = 1           # mm/pixel
#TODO: Add subroutine to make grey level intensity between 0 and 1

# Read and create aggregate object
data = reader.imageRead(fileName)
aggregate = Aggregate.Aggregate(fileName,data,pixelSize)

# Baseline data
boxSizeDEM, centres, radii = sdata.loadDEMboxsizeCentreRadius()
measure.measureBenchMarkSizeAndNormal(aggregate,radii,centres)

# Filter
print("\nIs the file in need of filtering(y/n)?:")
answer=input()
if answer=='y':
    filters.filterDenoiseNlm(aggregate)
    segment.greyLevelHistogram(aggregate)
else:
    aggregate.filteredGreyLevelMap = aggregate.greyLevelMap
    f = open("FilterDetails.txt","w+")
    f.write("No filtering parameters as no filter used")
    f.close()

# Binarize
segment.binarizeOtsu(aggregate)

# Clean-up
"""TODO: Closing holes in the binary maps. Needs to be iterative."""

#---Euclidean distance transform watershed
segment.euclidDistanceMap(aggregate)
segment.localMaximaMarkers(aggregate)
segment.topoWatershed(aggregate)

#---Particle size
measure.measureParticleSizeDistribution(aggregate)

#---Particle Contact 
measure.measureContactNormalsSpam(aggregate)

#---Figures


#%%---Writing log file


