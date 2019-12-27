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
import spam.datasets as sdata
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

# %% Objects

r = Reader.Reader()
f = Filter.Filter()
s = Segment.Segment()
m = Measure.Measure()
l = LemmeC.LemmeC()
w = Writer.Writer()
j = Jeeves.Jeeves()

# %% Reading data
calib = 11930/1000/1000 # mm/px

# Slice locations slices numbered 1 onwards: 
otcCenterSlice = round((20+1006)/2)
length = 5.5
upperSlice = otcCenterSlice + round( ( length / 2 ) / calib )
lowerSlice = otcCenterSlice - round( ( length / 2 ) / calib )

# Read entire data
file = r.readTiffSequence('/home/eg/github/PAC/pacInput/Deben1D/OTC/3D/0N',lowerSlice, upperSlice)

# %% Subregion extraction

# Center of sample rows/columns
csCenterRow = 450 # Y
csCenterCol = 506 # X
centerSlice = round( ( upperSlice - lowerSlice ) / 2 )


# %% 2 D50 cube subregion

twoD50 = 0.62*2
subUpperSlice = centerSlice + round( (twoD50 / 2) / calib )
subLowerSlice = centerSlice - round( (twoD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (twoD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (twoD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (twoD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (twoD50 / 2) / calib )
srTwoD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srTwoD50)

# Create Aggregates
two = Aggregate.Aggregate('two', srTwoD50, calib, 16)

# Filteration
f.filterDenoiseNlm(two)

# Binarization
s.binarizeOtsu( two )
s.fillholes( two )
s.removespecks( two )
two.binaryMap[ 0, :,: ] = 0
two.binaryMap[ two.binaryMap.shape[ 0 ] - 1, :, : ] = 0

# Euclid-Marker-Topo
s.euclidDistanceMap(two)
s.localhMaxima(two,4)
s.topoWatershed(two)
two.labelledMap = tiffy.imread('watershedSegmentation-edited.tif') # After manual editing
s.resetParticleList(two)

# Analysis
m.measureParticleSizeDistribution(two)

# %% 3 D50 cube subregion

threeD50 = 0.62*3
subUpperSlice = centerSlice + round( (threeD50 / 2) / calib )
subLowerSlice = centerSlice - round( (threeD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (threeD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (threeD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (threeD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (threeD50 / 2) / calib )
srThreeD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srThreeD50)

# Create Aggregates
three = Aggregate.Aggregate( 'three', srThreeD50, calib, 16 )

# Filteration
f.filterDenoiseNlm( three )

# Binarization
s.binarizeOtsu( three )
s.fillholes( three )
s.removespecks( three )
three.binaryMap[ 0, :,: ] = 0
three.binaryMap[ three.binaryMap.shape[ 0 ] - 1, :, : ] = 0

# Euclid-Marker-Topo
s.euclidDistanceMap(three)
s.localhMaxima(three,4)
s.topoWatershed(three)
three.labelledMap = tiffy.imread('watershedSegmentation-edited.tif') # After manual editing
s.resetParticleList(three)

# Analysis
m.measureParticleSizeDistribution(three)

# %% 4 D50 cube subregion

fourD50 = 0.62*4
subUpperSlice = centerSlice + round( (fourD50 / 2) / calib )
subLowerSlice = centerSlice - round( (fourD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (fourD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (fourD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (fourD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (fourD50 / 2) / calib )
srFourD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srFourD50)

# Create Aggregates
four = Aggregate.Aggregate( 'four', srFourD50, calib, 16 )

del srFourD50
gc.collect()

# Filteration
f.filterDenoiseNlm( four )

# Binarization
s.binarizeOtsu( four )
s.fillholes( four )
s.removespecks( four )
four.binaryMap[ 0, :,: ] = 0
four.binaryMap[ four.binaryMap.shape[ 0 ] - 1, :, : ] = 0

# Euclid-Marker-Topo
s.euclidDistanceMap( four )
s.localhMaxima( four, 4 )
s.topoWatershed( four )
four.labelledMap = tiffy.imread( 'watershedSegmentation-edited.tif' ) # After manual editing
s.resetParticleList( four )

# Analysis
m.measureParticleSizeDistribution( four )

# %% 5 D50 cube subregion

fiveD50 = 0.62*5
subUpperSlice = centerSlice + round( (fiveD50 / 2) / calib )
subLowerSlice = centerSlice - round( (fiveD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (fiveD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (fiveD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (fiveD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (fiveD50 / 2) / calib )
srFiveD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srFiveD50)

# Create Aggregates
five = Aggregate.Aggregate( 'five', srFiveD50, calib, 16 )

del srFiveD50
gc.collect()

# Filteration
f.filterDenoiseNlm( five )

# Binarization
s.binarizeOtsu( five )
s.fillholes( five )
s.removespecks( five )
five.binaryMap[ 0, :,: ] = 0
five.binaryMap[ five.binaryMap.shape[ 0 ] - 1, :, : ] = 0

# Euclid-Marker-Topo
s.euclidDistanceMap( five )
s.localhMaxima( five, 4 )
s.topoWatershed( five )
five.labelledMap = tiffy.imread( 'watershedSegmentation-edited.tif' ) # After manual editing
s.resetParticleList( five )

# Analysis
m.measureParticleSizeDistribution( five )





