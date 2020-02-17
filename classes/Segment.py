'''
Segment class
'''

# Importing libraries
from scipy.ndimage.morphology import distance_transform_edt as edt
from scipy.ndimage.morphology import binary_fill_holes as fillHoles
from scipy.ndimage.morphology import binary_opening as removeSpecks
from skimage.morphology import local_maxima as localMaxima
from skimage.morphology import h_maxima as hmax
from skimage.morphology import watershed as wsd
from skimage.filters import threshold_otsu
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt
import numpy as np


# %% Class

class Segment:

    def __init__( self ):
        print('\n-----------------')
        print('Segment activated')
        print('-----------------')

    def performEDTWS( self, filteredGLIMap, currentVoidRatio, outputFilesLocation, sampleName):
        '''
        Binarize
        EDM
        EDM peaks
        Topological watershed with EDM peaks
        '''

        print('Starting EDT-WS module')
        print('----------------------*\n')
        
        # Binarize
        print('Which binarization method to follow?')
        binarizationMethodToFollow = input('(1) Otsu threshold, (2) User input based threshold, ([3]) Density based threshold: ')
        
        # Otsu threshold
        if binarizationMethodToFollow == '1':
            binMap, gliThreshold = self.binarizeAccordingToOtsu( filteredGLIMap, outputFilesLocation, sampleName)
            voidRatioCT = self.calcVoidRatio( binMap )
        
        # User input based threshold
        elif binarizationMethodToFollow == '2':
            binMap, gliThreshold = self.binarizeAccordingToUserThreshold( filteredGLIMap, outputFilesLocation, sampleName )
            voidRatioCT = self.calcVoidRatio( binMap )
        
        # Density based threshold
        else:
            binMap, gliThreshold = self.binarizeAccordingToDensity(filteredGLIMap, currentVoidRatio, outputFilesLocation, sampleName )
            voidRatioCT = self.calcVoidRatio( binMap )

        # EDM
        edMap = self.obtainEuclidDistanceMap( binMap )
         
        # Markers
        edPeakMrkrMap = self.obtainLocalMaximaMarkers( edMap )
        
        # Labelled Map
        labelledMap = self.obtainLabelledMapUsingWaterShedAlgorithm( binMap, edMap, edPeakMrkrMap )
        
        # Correction of labelled map
        lblCorrectionMethod, correctedLabelledMap = self.fixErrorsInSegmentation( labelledMap )

        # Returns
        return gliThreshold, binMap, voidRatioCT, edMap, edPeakMrkrMap, lblCorrectionMethod, correctedLabelledMap

    def binarizeAccordingToOtsu( self, gliMapToBinarize, outputLocation, sampleName ):       
        print('\nRunning Otsu Binarization')
        print('----------------------------*')
        otsuThreshold = threshold_otsu( gliMapToBinarize )      
        binaryMap = np.zeros_like( gliMapToBinarize )        
        
        binaryMap[ np.where( gliMapToBinarize > otsuThreshold ) ] = 1
        binaryMap = binaryMap.astype( int )
        e1 = self.calcVoidRatio( binaryMap )
        
        binaryMap = self.fillHoles( binaryMap )
        e2 = self.calcVoidRatio( binaryMap )        
        
        binaryMap = self.removeSpecks( binaryMap )      
        e3 = self.calcVoidRatio( binaryMap ) 

        print( 'Global Otsu threshold = %f' % otsuThreshold )
        print( 'Void ratio after threshold = %f' % e1 )  
        print( 'Void ratio after filling holes = %f' % e2 )   
        print( 'Void ratio after removing specks = %f' % e3 )       
        
        # Saving files
        otsuThresholdFileName = outputLocation + sampleName + '-otsuThresholdDetails.txt'
        f = open( otsuThresholdFileName, 'w+' )        
        f.write( 'Global Otsu threshold = %f' % otsuThreshold )
        f.write( '\nVoid ratio after threshold = %f' % e1 )  
        f.write( '\nVoid ratio after filling holes = %f' % e2 )   
        f.write( '\nVoid ratio after removing specks = %f' % e3 )
        f.close()

        print('Output file saved as ' + otsuThresholdFileName)
        print('---------------------*')

        return binaryMap, otsuThreshold 
              
    def binarizeAccordingToUserThreshold( self, gliMapToBinarize, outputLocation, sampleName):        
        userThreshold = int( input( 'Enter user threshold: ' ) )        
        binaryMap = np.zeros_like( gliMapToBinarize )        
        binaryMap[ np.where( gliMapToBinarize > userThreshold ) ] = 1
        binaryMap = binaryMap.astype( int )
        e1 = self.calcVoidRatio( binaryMap ) 
                
        # Filling Holes
        binaryMap = fillHoles( binaryMap )      
        e2 = self.calcVoidRatio( binaryMap )    
  
        # Removing specks
        binaryMap = self.removeSpecks( binaryMap )
        e3 = self.calcVoidRatio( binaryMap )     

        print( 'Global User threshold = %f' % userThreshold )
        print( 'Void ratio after threshold = %f' % e1 )  
        print( 'Void ratio after filling holes = %f' % e2 )   
        print( 'Void ratio after removing specks = %f' % e3 )    
                
        # Saving files
        userThresholdFileName = outputLocation + sampleName + '-userThresholdDetails.txt'
        f = open( userThresholdFileName, 'w+' )        
        f.write( 'Global user threshold = %f' % userThreshold )
        f.write( '\nVoid ratio after threshold = %f' % e1 )  
        f.write( '\nVoid ratio after filling holes = %f' % e2 )   
        f.write( '\nVoid ratio after removing specks = %f' % e3 )
        f.close() 

        print('Output file saved as ' + userThresholdFileName)
        print('---------------------*')

        return binaryMap, userThreshold

    def binarizeAccordingToDensity( self, gliMapToBinarize, knownVoidRatio, outputLocation, sampleName):      
        otsuBinMap, otsuThreshold = self.binarizeAccordingToOtsu(gliMapToBinarize, outputLocation, sampleName)
        currentThreshold = otsuThreshold
        currentVoidRatio = self.calcVoidRatio( otsuBinMap )
        targetVoidRatio = knownVoidRatio
        greyLvlMap = gliMapToBinarize
            
        tolerance = 0.0001
        deltaVoidRatio = targetVoidRatio - currentVoidRatio
            
        absDeltaThreshold = 0
        maxAbsDeltaThreshold = 500
            
        iterationNum = 1
        maxIterations = 50
            
        userThresholdStepsTextFileName = outputLocation +  sampleName + '-userThresholdDensitySteps.txt'
        f = open( userThresholdStepsTextFileName, "w+")
        f.write( "Steps of density based user threshold\n\n" )
        f.close()

        while( abs( deltaVoidRatio ) > tolerance and iterationNum <= maxIterations ):
            incrementSign = deltaVoidRatio/abs(deltaVoidRatio)  
               
            if abs( deltaVoidRatio ) > 0.100:
                absDeltaThreshold = 500
                
            elif abs( deltaVoidRatio ) > 0.05:
                absDeltaThreshold = 250
                
            elif abs( deltaVoidRatio ) > 0.01:
                absDeltaThreshold = 100
                
            elif abs( deltaVoidRatio ) > 0.005:
                absDeltaThreshold = 50
                
            elif abs( deltaVoidRatio ) > 0.001:
                absDeltaThreshold = 25
                
            elif abs( deltaVoidRatio ) > 0.0005:
                absDeltaThreshold = 10
                
            else:
                absDeltaThreshold = 1
                
            currentThreshold = currentThreshold + absDeltaThreshold*incrementSign
            currentBinaryMap = np.zeros_like( greyLvlMap )        
            currentBinaryMap[ np.where( greyLvlMap > currentThreshold ) ] = 1
                
            currentBinaryMap = self.fillHoles( currentBinaryMap )
            currentBinaryMap = self.removeSpecks( currentBinaryMap )
                
            currentVoidRatio = self.calcVoidRatio( currentBinaryMap )
            deltaVoidRatio = targetVoidRatio - currentVoidRatio
                
            print( '\nIteration %d:' % iterationNum )
            print( 'Current threshold: %d' % round( currentThreshold ) )
            print( 'Otsu Threshold: %d' % round( otsuThreshold ) )
            print( 'Current void ratio: %0.5f' % round( currentVoidRatio,5 ) )
            print( 'Target void ratio: %0.5f' % round( targetVoidRatio,5 ) )
            print( 'Target - Current void ratio: %0.5f' % round( deltaVoidRatio, 5 ) )

            f = open( userThresholdStepsTextFileName, 'a' )
            f.write( "\nIteration %d:\n" % iterationNum)
            f.write( "Current threshold: %d\n" % round( currentThreshold ) )
            f.write( "Otsu Threshold: %d\n" %  round( otsuThreshold ) )
            f.write( "Current void ratio:  %0.5f\n" % round( currentVoidRatio,5 ) )
            f.write( "Target void ratio: %0.5f\n" % round( targetVoidRatio,5 ) )
            f.write( "Target - Current void ratio: %0.5f\n" % round( deltaVoidRatio, 5 ) )
            f.close()

            iterationNum = iterationNum + 1
            
            
        densityBasedThresholdTextFileName = outputLocation +  sampleName + '-userThresholdDensityBased.txt'
        f = open( densityBasedThresholdTextFileName, "w+" )               
        f.write( "Global density based threshold = %f\n" % currentThreshold )
        f.close()

        return currentBinaryMap, currentThreshold
              
    def calcVoidRatio(self, binaryMapforVoidRatioCalc):
        '''
        Parameters
        ----------
        binaryMap : numpy array with 1 and 0
        
        Returns
        -------
        void ratio of the binary map assumin 1 is particle
        '''
        volSolids = binaryMapforVoidRatioCalc.sum()
        lenZ = binaryMapforVoidRatioCalc.shape[ 0 ]
        lenY = binaryMapforVoidRatioCalc.shape[ 1 ]
        lenX = binaryMapforVoidRatioCalc.shape[ 2 ]
        volTotal =  lenZ * lenY * lenX 
        currentVoidRatio = ( volTotal - volSolids ) / volSolids
        return currentVoidRatio
        
    def fillHoles( self, oldBinaryMap ):
        '''
        
        Parameters
        ----------
        oldBinaryMap : Numpy array of binary data
            obtained from binarization

        Returns
        -------
        newBinaryMap holes in the map is filled

        '''
        newBinaryMap = fillHoles( oldBinaryMap )
        return newBinaryMap.astype(int) 
                  
    def removeSpecks( self, oldBinaryMap ):
        '''
        
        Parameters
        ----------
        oldBinaryMap : Numpy array of binary data
            obtained from binarization

        Returns
        -------
        newBinaryMap with white specks in the map removed

        '''
        newBinaryMap = removeSpecks( oldBinaryMap )   
        return newBinaryMap.astype(int)
       
    def obtainEuclidDistanceMap( self, binaryMapForEDM ):      
        print('\nFinding Euclidian distance map (EDM)')
        print('------------------------------------*')
        edMap = edt( binaryMapForEDM )             
        print( "EDM Created" )        
        return edMap

    def obtainLocalMaximaMarkers( self, edMapForPeaks ):    
        '''
        Fix this to be better at obtaining markers
        Expected number of particles?
        Particle size?
        '''
        print('\nObtaining peaks of EDM...')
        print('------------------------------*')  
        h = int( input( 'Enter the minimum height for a peak (px): ') )
        print( 'Finding local maxima in EDM' )
                
        edmPeakMarkers = hmax( edMapForPeaks, h ).astype(int)    
        
        print( 'Found local maximas in EDM (greater than %i px)' % h )

        print( 'Resetting count of peaks' )
        count=0    
        for frame in range( 0, edmPeakMarkers.shape[ 0 ] ):            
            for row in range( 0, edmPeakMarkers.shape[ 1 ] ):                
                for col in range( 0, edmPeakMarkers.shape[ 2 ] ):                    
                    if edmPeakMarkers[ frame ][ row ][ col ] == 1:                        
                        edmPeakMarkers[ frame ][ row ][ col ] = count + 1                        
                        count = count + 1
            print( 'Processed ' + str(frame + 1) + ' out of ' + str(edmPeakMarkers.shape[ 0 ]) + ' slices' )
        
        return edmPeakMarkers

    def obtainLabelledMapUsingWaterShedAlgorithm(self, binaryMap, euclidDistMap, edPeaksMap ):
        '''
        Returns labelled map
        '''
        print( '\nSegmenting particles by topological watershed' ) 
        binMask = binaryMap.astype( bool )
        invertedEdMap = -euclidDistMap
        labelledMap = wsd( invertedEdMap, markers = edPeaksMap, mask = binMask )
        print( 'Watershed segmentation complete' )
        
        return labelledMap

    def fixErrorsInSegmentation(self, aggregate):
        '''
        Parameters
        ----------
        aggregate : Object
            

        Returns
        -------
        Updates aggregate object to correct over-segmentation and fixes the particle lists

        '''
        checkCorrectionMethod = input( 'Manual or automated error correction? ["m"/a]: ' )

        if checkCorrectionMethod.lower() == 'a':
            print( 'Not coded yet' )

        else:
            print( '-------------------------------------*' )
            print( '\n\nEntering manual correction method' )
            print( 'Carry out manual correction in ImageJ' )
            print( 'Once manual correction is complete, save corrected labelled file as sandName-updatedWatershedSegmentation.tiff in output folder' )
            correctedWatershedFileName = input( 'Enter name of corrected file: ')
            completeCorrectedWatershedFileName = aggregate.dataOutputDirectory + correctedWatershedFileName
            aggregate.labelledMap = tiffy.imread(completeCorrectedWatershedFileName)

        return labelCorrectionMethod, correctedLabelMap
        
    def resetParticleList(self, aggregate):
        '''
        Parameters
        ----------
        aggregate : Object
            

        Returns
        -------
        Updates particle list after manual correction of segmentation
 
        '''
        aggregate.numberOfParticles = int( aggregate.labelledMap.max() )        
        aggregate.particleList = []         
        aggregate.particleList.append(None) 
        
        for i in range( 1, aggregate.numberOfParticles + 1 ):          
            numberOfParticleVoxel = np.where( aggregate.labelledMap == i)[ 0 ].shape[ 0 ]            
            locData = np.zeros( ( numberOfParticleVoxel, 3 ) )                            
            
            for j in range(0, 3):              
                locData[ :, j ] = np.where( aggregate.labelledMap == i )[ j ]                
            
            p = Particle.Particle( i, numberOfParticleVoxel, locData )                             
            aggregate.particleList.append( p )

