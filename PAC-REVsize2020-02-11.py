# -*- coding: utf-8 -*-

"""
    This is for analysis of 1D compression data scans at no load to obtain 
        1. REV size
        2. Particle shape parameter
"""

#%% Imports

print("\n\nStarting Particle analysis code (PAC) for REV size analysis------------*\n\n")

# General imports
import classes.Aggregate as Aggregate

# Verb imports
import classes.Reader as Reader
import classes.Filter as Filter
import classes.Segment as Segment
import classes.Measure as Measure
import classes.LemmeC as LemmeC
import classes.Jeeves as Jeeves

# Objects
r = Reader.Reader()
f = Filter.Filter()
s = Segment.Segment()
m = Measure.Measure()
l = LemmeC.LemmeC()
j = Jeeves.Jeeves()

# %% SuperCube analysis

j.checkFolderLocations()
j.checkREVSampleDetails()
gliStack = r.readTiffSequence( j.tiffFileLocation, j.superCubeCenterSlice, j.superCubeCenterRow, j.superCubeCenterCol, j.superCubeEdgeLengthforREVAnalysis, j.superCubeCalib )
superCube = Aggregate.Aggregate( j.sandName, gliStack, j.superCubeCalib, 16, j.voidRatio, True )
f.checkForFiltration(superCube)    
s.binarizeOtsu( superCube )
s.resetBinarizationAccordingToDensity( superCube )


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

