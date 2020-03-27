'''

'''

import numpy as np
import skimage.external.tifffile as tiffy

from pac import Particle as Particle
from pac import Reader as Reader

class Aggregate:

    def __init__( self, sampleName, bitDpth, sampleCalib, cubeEdgeLength, centerZ, centerY, centerX, d50Known, voidRatioKnown, dataFromTiffStack, tiffFileLocation, unfilteredGLIData):
        self.sampleLocation = tiffFileLocation
        self.sampleName = sampleName
        self.bitDepth = bitDpth
        self.gliMax = 2 ** bitDpth - 1
        self.calib = sampleCalib

        self.sampleEdgeLength = cubeEdgeLength
        self.sampleCenterZ = centerZ
        self.sampleCenterY = centerY
        self.sampleCenterX = centerX

        self.originalD50FromSieve = d50Known
        self.currentD50FromSieve = 0.0

        self.currentvoidRatioMeasured = voidRatioKnown
        self.currentvoidRatioFromCT = 0.0

        self.dataTakenFromTiffStack = dataFromTiffStack
        self.locationOfDataFile = tiffFileLocation

        self.greyLevelMap = unfilteredGLIData
        self.filteredGreyLevelMap = self.greyLevelMap
        self.imageNoise = np.zeros_like( self.filteredGreyLevelMap )
        self.greyLevelHistogram = np.zeros( ( self.gliMax, 2 ) )
        self.filteredGreyLevelHistogram = np.zeros( ( self.gliMax, 2 ) )

        self.gliThreshold = 0
        self.binaryMap = np.zeros_like( self.greyLevelMap )

        self.euclidDistanceMap = np.zeros_like( self.greyLevelMap )
        self.edPeakMarkers = np.zeros_like( self.greyLevelMap )

        self.completeSegmentationMethod = ''
        self.labelledMap = np.zeros_like( self.greyLevelMap )
        self.particleList = [ ]
        self.numberOfParticles = 0

        self.particleSizeDataSummaryFromEDTWS= np.zeros( ( 5000, 6 ) ) 
        self.grainSizeDistributionEquivalentSphereFromEDTWS = np.zeros( ( 5000, 2 ) )
        self.grainSizeDistributionFeretMaxFromEDTWS = np.zeros( ( 5000, 2 ) )
        self.grainSizeDistributionFeretMinFromEDTWS = np.zeros( ( 5000, 2 ) )
        self.grainSizeDistributionFeretMedFromEDTWS = np.zeros( ( 5000, 2 ) )


        print("\nAggregate activated")
        print('--------------------*')

