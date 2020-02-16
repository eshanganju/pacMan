'''

'''

import numpy as np
import skimage.external.tifffile as tiffy

from classes import Particle as Particle
from classes import Reader as Reader

class Aggregate:

    def __init__( sampleName, sampleCalib, cubeEdgeLength, samplCenterZ, sampleCenterY, sampleCenterZ, d50Known, voidRatioKnown, tiffFileFolderLocation):
        self.dataOutputDirectory = dataStorageDirectory
        self.sampleLocation = sampleLocation
        self.sampleName = sampleName
        self.pixelSize = pxlSize
        self.bitDepth = bitDpth
        self.GLIMax = 2 ** bitDpth - 1
        self.calib = mmToPxcalibration
        
        
        self.greyLevelMap = greyLvlMap
        self.filteredGreyLevelMap = greyLvlMap
        self.imageNoise = np.zeros_like( self.filteredGreyLevelMap )
        self.binaryMap = np.zeros_like( self.greyLevelMap )
        self.binaryMapUser = np.zeros_like( self.greyLevelMap )
        self.euclidDistanceMap = np.zeros_like( self.greyLevelMap )
        self.markers = np.zeros_like( self.greyLevelMap )
        self.labelledMapFromEDTWS = np.zeros_like( self.greyLevelMap )      
        self.particleList = [ ]         
        self.particleList.append( None ) # Background taken as 0th particle

        self.superCubeCenterSlice = int( 0 )
        
        self.globalOtsuBasedThreshold = 0
        self.globalUserBasedThreshold = 0
        self.globalDensityBasedThreshold = 0

        self.numberOfParticles = 0      
        
        self.particleSizeDataSummary = np.zeros( ( 5000, 6 ) ) 
        self.grainSizeDistributionEquivalentSphere = np.zeros( ( 5000, 2 ) )     
        self.grainSizeDistributionFeretMax = np.zeros( ( 5000, 2 ) )         
        self.grainSizeDistributionFeretMin = np.zeros( ( 5000, 2 ) )             
        self.grainSizeDistributionFeretMed = np.zeros( ( 5000, 2 ) )             
        
        self.greyLevelHistogram = np.zeros( ( 1000, 2 ) )                    
        self.filteredGreyLevelHistogram = np.zeros( ( 1000, 2 ) )
               
        self.measuredD50 = 0
        self.ctD50 = 0

        self.measuredVoidRatio = 0       
        self.ctVoidRatio = 0
        
        createGLIFile = input('Do you want to save a copy of the GLI[y/"n"]')
        if createGLIFile.lower() == 'y':
            unfilteredImageName = str(self.dataOutputDirectory) + str( self.fileName ) + '.tiff'
            tiffy.imsave( unfilteredImageName, self.greyLevelMap )
            print( "\nUnfiltered GLI file saved as " + imageName )         
        else:
            print('File not saved...')
        
        print("\nAggregate activated")
        print('--------------------*')

