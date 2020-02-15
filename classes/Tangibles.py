'''

'''

import numpy as np
import skimage.external.tifffile as tiffy

class Aggregate:

    def __init__( self, dataStorageDirectory, sampleLocation, sampleName, pxlSize, bitDpth, mmToPxcalibration, greyLvlMap):
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


class Particle:

    def __init__(self, num, numVox, locData):
        self.index = num
        self.volume = numVox
        self.locationData = locData
        self.covarianceMatrix = np.zeros((3,3))
        self.eigenValue = np.zeros(3)
        self.eigenVector = np.zeros((3,3))
        self.meanZ = 0
        self.meanY = 0
        self.meanX = 0
        self.centeredLocationData = np.zeros_like(self.locationData)
        self.rotatedCenteredLocationData = np.zeros_like(self.locationData)
        self.feretMax = 0
        self.feretMin = 0
        self.feretMed = 0
        self.equivalentSphereDiameter = 0
        print("Particle No.",self.index,"created")

