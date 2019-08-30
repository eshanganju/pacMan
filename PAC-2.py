# -*- coding: utf-8 -*-

"""
This is the second version of the Particle analysis code

The following are attempted with trepidation: 
  1. Reading tiff files containing CT scan data in to np arrays
  2. Filtering the CT scan data to reduce noise
  3. Binarization and segmentation of the data to separate out individual particles
  4. Computing particle size and morphology distribution
  5. Location interparticle contact
  6. Determining contact normal distribution

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

# General
import numpy as np
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt

# Made up
import Reader
import Particle
import Aggregate
import Filter
import Segment
import Measure
import LemmeC
import Writer

# Activating objects
reader = Reader.Reader()
filters = Filter.Filter()
segment = Segment.Segment()
measure = Measure.Measure()
lemmec = LemmeC.LemmeC()
writer = Writer.Writer()

#---Filename
fileName = 'Box4A.tiff'
pixelSize = 1
f = open("logFile.txt","w+")
f.write("File name: %s " % fileName)


#---Read
data = reader.imageRead(fileName)
aggregate = Aggregate.Aggregate(data,pixelSize)


#---Filter
patchSize,patchDist,cutOff = filters.filterDenoiseNlm(aggregate)
segment.greyLevelHistogram(aggregate)
f.write("\n\nFilter parameters---------------------*\n")
f.write("Patch size = %f\n" % patchSize)
f.write("Patch distance = %f\n" % patchDist)
f.write("Cut-off intensity = %f\n" % cutOff)

#---Binarize
segment.binarizeOtsu(aggregate)

#---Clean-up


#---ED map
segment.euclidDistanceMap(aggregate)

#---Markers
segment.localMaximaMarkers(aggregate)

#---Segment
segment.topoWatershed(aggregate)

#---Particle size
measure.measureParticleSizeDistribution(aggregate)

#---Particle Contact
measure.measureContactNormalsSpam(aggregate)


#---Figures


f.close()

