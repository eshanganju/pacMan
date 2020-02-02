# -*- coding: utf-8 -*-

"""
    This is for analysis of 1D compression data scans at no load to obtain 
        1. REV size
        2. Particle shape parameter
"""

# %% Imports

print("\n\nStarting Particle analysis code (PAC)------------*")

# General imports
import numpy as np
import skimage.external.tifffile as tiffy
from skimage import io
import matplotlib.pyplot as plt
import gc

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

# Objects

r = Reader.Reader()
f = Filter.Filter()
s = Segment.Segment()
m = Measure.Measure()
l = LemmeC.LemmeC()
w = Writer.Writer()
j = Jeeves.Jeeves()

# Subregion extraction
'''
    OTC:0.72
    OGF: 0.62
    2QR: 0.71
'''
D50 = 0.72

# Reading data
calib = 11930/1000/1000 # mm/px

#  Center slice location: 
otcCenterSlice = round((20+1006)/2)

# Center of sample rows/columns
csCenterRow = 450 # Y
csCenterCol = 506 # X

# Center slice
length = 5.5
upperSlice = otcCenterSlice + round( ( length / 2 ) / calib )
lowerSlice = otcCenterSlice - round( ( length / 2 ) / calib )
upperRow = csCenterRow + round( (length / 2) / calib )
lowerRow = csCenterRow - round( (length / 2) / calib )
upperCol = csCenterCol + round( (length / 2) / calib )
lowerCol = csCenterCol - round( (length / 2) / calib )

# Read entire data
file = r.readTiffSequence('C:/Users/eganj/gitHub/pac/OTC-0N',lowerSlice, upperSlice-1)[:,lowerRow:upperRow,lowerCol:upperCol]

# %% Determination of user threshold for binarization (using a 5.5 mm cube volume)

'''
    Find out the GLI threshold that results in same void ratio as measured value from sample
    The volume of cube used is 5.5 mm
    Void ratios measured from test: 
        OTC: 0.5404
        OGF: 0.6348
        2QR: 0.7337
'''

# Create Aggregates 
superCube = Aggregate.Aggregate('5.5mm cube', file, calib, 16)
tiffy.imsave('superCube.tiff', superCube.greyLevelMap)

# Filteration
f.filterDenoiseNlm(superCube)
r.plotGLI(superCube.filteredGreyLevelMap)

superCube.filteredGreyLevelMap = tiffy.imread('C:/Users/eganj/Google Drive/(02) Research/(07) EIDPS/(02) Data/(05) Tomo/(02) Physics - nCT/(05) 1D compression study/(02) REV analysis/density based binarization/OTC/superCube/Filtered.tiff')

# Binarization - OTSU
s.binarizeOtsu( superCube )
print('Otsu Threshold = %f' % superCube.globalOtsuThreshold)

# Refine OTSU
thresholdUser = 20256
s.resetOtsuBinarizationAccordingToUser(superCube,thresholdUser)

# %% Cube subregion

# Cubical volume
nD50 = 3
sublength = nD50*D50

superCubeCenterSlice = (superCube.greyLevelMap.shape[0])//2 
superCubeCenterRow = (superCube.greyLevelMap.shape[1])//2  # Y
superCubeCenterCol = (superCube.greyLevelMap.shape[2])//2  # Z

subUpperSlice = superCubeCenterSlice + round( (sublength / 2) / calib )
subLowerSlice = superCubeCenterSlice - round( (sublength / 2) / calib ) 
subUpperRow = superCubeCenterRow + round( (sublength / 2) / calib )
subLowerRow = superCubeCenterRow - round( (sublength / 2) / calib )
subUpperCol = superCubeCenterCol + round( (sublength / 2) / calib )
subLowerCol = superCubeCenterCol - round( (sublength / 2) / calib )

gli = superCube.greyLevelMap[ subLowerSlice : subUpperSlice , subLowerRow : subUpperRow , subLowerCol : subUpperCol ]
fgli = superCube.filteredGreyLevelMap[ subLowerSlice : subUpperSlice , subLowerRow : subUpperRow , subLowerCol : subUpperCol ]

r.plotGLI(gli)
r.plotGLI(fgli)

# Create Aggregate
cube = Aggregate.Aggregate( nD50, gli, calib, 16 )
cube.filteredGreyLevelMap = fgli

# Binarize Otsu
s.binarizeOtsu( cube )
print('Otsu Threshold = %f' % cube.globalOtsuThreshold)

    
# Refine OTSU
s.resetOtsuBinarizationAccordingToUser(cube,thresholdUser)

# Euclid-Marker-Topo
s.euclidDistanceMap( cube )

# Particle markers
s.localhMaxima( cube, 4 )

# Watershed
s.topoWatershed( cube )

# Update particles
cube.labelledMap = tiffy.imread( 'watershedSegmentation-edited3-NE.tif' ) # After manual editing
s.resetParticleList( cube )

# Analysis
m.measureParticleSizeDistribution( cube )

