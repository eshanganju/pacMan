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
file = r.readTiffSequence('/home/eg/github/PAC/pacInput/Deben1D/3D/0N',lowerSlice, upperSlice)

# %% Subregion extraction

# Center of sample rows/columns
csCenterRow = 450 # Y
csCenterCol = 506 # X
centerSlice = round( ( upperSlice - lowerSlice ) / 2 )

# 2 D50 cube subregion
twoD50 = 0.62*2
subUpperSlice = centerSlice + round( (twoD50 / 2) / calib )
subLowerSlice = centerSlice - round( (twoD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (twoD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (twoD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (twoD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (twoD50 / 2) / calib )
srTwoD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srTwoD50)

#%% 
# 3 D50 cube subregion
threeD50 = 0.62*3
subUpperSlice = centerSlice + round( (threeD50 / 2) / calib )
subLowerSlice = centerSlice - round( (threeD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (threeD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (threeD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (threeD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (threeD50 / 2) / calib )
srThreeD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srThreeD50)

# 4 D50 cube subregion
fourD50 = 0.62*4
subUpperSlice = centerSlice + round( (fourD50 / 2) / calib )
subLowerSlice = centerSlice - round( (fourD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (fourD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (fourD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (fourD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (fourD50 / 2) / calib )
srFourD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srFourD50)

# 5 D50 cube subregion
fiveD50 = 0.62*5
subUpperSlice = centerSlice + round( (fiveD50 / 2) / calib )
subLowerSlice = centerSlice - round( (fiveD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (fiveD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (fiveD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (fiveD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (fiveD50 / 2) / calib )
srFiveD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srFiveD50)

# 6 D50 cube subregion
sixD50 = 0.62*6
subUpperSlice = centerSlice + round( (sixD50 / 2) / calib )
subLowerSlice = centerSlice - round( (sixD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (sixD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (sixD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (sixD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (sixD50 / 2) / calib )
srSixD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srSixD50)

# 7 D50 cube subregion
sevenD50 = 0.62*7
subUpperSlice = centerSlice + round( (sevenD50 / 2) / calib )
subLowerSlice = centerSlice - round( (sevenD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (sevenD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (sevenD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (sevenD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (sevenD50 / 2) / calib )
srSevenD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srSevenD50)

# 8 D50 cube subregion
eightD50 = 0.62*8
subUpperSlice = centerSlice + round( (eightD50 / 2) / calib )
subLowerSlice = centerSlice - round( (eightD50 / 2) / calib ) 
csSubUpperRow = csCenterRow + round( (eightD50 / 2) / calib )
csSubLowerRow = csCenterRow - round( (eightD50 / 2) / calib )
csSubUpperCol = csCenterCol + round( (eightD50 / 2) / calib )
csSubLowerCol = csCenterCol - round( (eightD50 / 2) / calib )
srEightD50 = file[ subLowerSlice : subUpperSlice , csSubLowerRow : csSubUpperRow , csSubLowerCol : csSubUpperCol ]
r.plotGLI(srEightD50)

del file

# %% Segmentation

# Create Aggregates
two = Aggregate.Aggregate('two', srTwoD50, calib, 16)

# Filteration
f.filterDenoiseNlm(two)

# Binarization
s.binarizeOtsu(two)
s.fillholes(two)
s.removespecks(two)

# Euclid-Marker-Topo
s.euclidDistanceMap(two)
s.localhMaxima(two,12)
s.topoWatershed(two)

#%% Purg

s.euclidDistanceMap(two)
