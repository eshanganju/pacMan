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

#---General Classes
import numpy as np
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt
import spam.datasets as sdata

#---Unique Classes
import Reader
import Aggregate
import Filter
import Segment
import Measure
import LemmeC
import Writer
import Jeeves

