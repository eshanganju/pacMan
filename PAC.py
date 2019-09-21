# -*- coding: utf-8 -*-

"""
This is the second version of the Particle analysis code

The following are attempted with trepidation: 
  1. Read tiff files containing CT scan data in to np arrays and normalize grey scale to [0,1]
  2. Filter CT scan data to reduce noise using automated NLM loops
  3. Segment CT data using:  
    3a. Traditional topographical watershed segmentation
    3b. Power segmentation
    3c. Level-set segmentation
  4. Compute particle size and morphology distribution
  5. Locate interparticle contact
  6. Determine contact normal distribution

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

#---General Classes
import numpy as np
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt
import spam.datasets as sdata

#---Unique Classes
import Reader
import Particle
import Aggregate
import Filter
import Segment
import Measure
import LemmeC
import Writer

#---Objects

print("\n\nStarting Particle analysis code (PAC)------------*")
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

#---Baseline data
boxSizeDEM, centres, radii = sdata.loadDEMboxsizeCentreRadius()
measure.measureBenchMarkSizeAndNormal(aggregate,radii,centres)

#---Filter
patchSize,patchDist,cutOff = filters.filterDenoiseNlm(aggregate)
f.write("\n\nFilter parameters---------------------*\n") 
f.write("Patch size = %f\n" % patchSize) 
f.write("Patch distance = %f\n" % patchDist)
f.write("Cut-off intensity = %f\n" % cutOff)

#---Binarize
segment.binarizeOtsu(aggregate)

#---Clean-up
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


f.close()

