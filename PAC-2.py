# -*- coding: utf-8 -*-

"""
This is the second version of the Particle analysis code

The following are attempted with trepidation: 
  1. Reading tiff files containing CT scan data in to np arrays
  2. Filtering the CT scan data to reduce noise
  3. Segmenting the data to separate out individual particles
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

fileName = 'Box4A.tiff'
pixelSize = 1

data = reader.imageRead(fileName)
aggregate = Aggregate.Aggregate(data,pixelSize)

###
plt.figure()
cut=aggregate.greyLevelMap.shape[0]//2
plt.imshow(aggregate.greyLevelMap[cut])
plt.show()
###


