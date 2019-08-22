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
    3. Filter
    4. Segment
    5. Measure
    6. Visualize
"""

# importing important stuff

import numpy as np

