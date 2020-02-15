'''

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

from classes import Tangibles as Tangibles

# %% Class

class Segment:

    def __init__( self ):
        print( "Segmenter activated" )


    def binarizeAccordingToOtsu( self, aggregate ):       
        print('\nRunning Otsu Binarization')
        greyLvlMap = aggregate.filteredGreyLevelMap
        aggregate.globalOtsuThreshold = threshold_otsu( greyLvlMap )      
        aggregate.binaryMap = np.zeros_like( greyLvlMap )        
        aggregate.binaryMap[ np.where( greyLvlMap > aggregate.globalOtsuThreshold ) ] = 1
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e1 = self.calcVoidRatio(aggregate.binaryMap)
        print( 'Void ratio after Otsu binarization= %f' % e1 )
          
        # Filling Holes
        aggregate.binaryMap = fillHoles(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e2 = self.calcVoidRatio(aggregate.binaryMap)        
        print('Void ratio after filling holes = %f' % e2)
        
        # Removing specks
        aggregate.binaryMap = removeSpecks(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )        
        e3 = self.calcVoidRatio(aggregate.binaryMap)        
        aggregate.ctVoidRatio = e3       
        print('Void ratio after removing specks = %f' % e3)
        
        # Saving files
        print( "Completed binarization using Otsu threshold of %d" % aggregate.globalOtsuThreshold)                
        otsuBinaryFileName = aggregate.dataOutputDirectory + aggregate.fileName+'-otsuBinary.tiff'
        tiffy.imsave( otsuBinaryFileName, aggregate.binaryMap )        
        otsuTextFileName = aggregate.dataOutputDirectory + aggregate.fileName + '-OtsuThreshold.txt'
        f = open( otsuTextFileName, "w+" )        
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalOtsuThreshold )        
               
    def binarizeAccordingToUserThreshold( self, aggregate, userThreshold):        
        print('\nRunning Otsu Binarization')
        greyLvlMap = aggregate.filteredGreyLevelMap
        aggregate.globalOtsuThreshold = threshold_otsu( greyLvlMap )      
        aggregate.binaryMap = np.zeros_like( greyLvlMap )        
        aggregate.binaryMap[ np.where( greyLvlMap > aggregate.globalOtsuThreshold ) ] = 1
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e1 = self.calcVoidRatio(aggregate.binaryMap)
        print( 'Void ratio after Otsu binarization= %f' % e1 )
               
        # Filling Holes
        aggregate.binaryMap = fillHoles(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e2 = self.calcVoidRatio(aggregate.binaryMap)        
        print('Void ratio after filling holes = %f' % e2)
        
        # Removing specks
        aggregate.binaryMap = removeSpecks(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )        
        e3 = self.calcVoidRatio(aggregate.binaryMap)        
        aggregate.ctVoidRatio = e3       
        print('Void ratio after removing specks = %f' % e3)
        
        # Saving files
        print( "Completed binarization using Otsu threshold of %d" % aggregate.globalOtsuThreshold)                
        otsuBinaryFileName = aggregate.dataOutputDirectory + aggregate.fileName+'-otsuBinary.tiff'
        tiffy.imsave( otsuBinaryFileName, aggregate.binaryMap )        
        otsuTextFileName = aggregate.dataOutputDirectory + aggregate.fileName + '-OtsuThreshold.txt'
        f = open( otsuTextFileName, "w+" )        
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalOtsuThreshold )

        # Resetting binarization threshold    
        aggregate.globalUserThreshold = userThreshold        
        aggregate.binaryMap = np.zeros_like( greyLvlMap )        
        aggregate.binaryMap[ np.where( greyLvlMap > aggregate.globalUserThreshold ) ] = 1
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e4 = self.calcVoidRatio(aggregate.binaryMap) 
        print( 'void ratio after binarization= %f' % e4 )
                
        # Filling Holes
        aggregate.binaryMap = fillHoles(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )       
        e5 = self.calcVoidRatio(aggregate.binaryMap)    
        print('void ratio after filling holes = %f' % e5)
                
        # Removing specks
        aggregate.binaryMap = removeSpecks(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e6 = self.calcVoidRatio(aggregate.binaryMap)     
        aggregate.ctVoidRatio = e6      
        print('void ratio after removing specks = %f' % e)
                
        # Saving files
        print( "Completed binarization using user input" )  
        userBinaryFileName = aggregate.dataOutputDirectory + aggregate.fileName +'-userBinary.tiff'              
        tiffy.imsave( userBinaryFileName, aggregate.binaryMapUser )        
        
        userThresholdTextFileName = aggregate.dataOutputDirectory + aggregate.fileName + '-userThreshold.txt'
        f = open( userThresholdTextFileName, "w+" )               
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalUserThreshold )        
        f.close()        

    def binarizeAccordingToDensity( self, aggregate ):
        print('\nRunning Otsu Binarization')
        greyLvlMap = aggregate.filteredGreyLevelMap
        aggregate.globalOtsuThreshold = threshold_otsu( greyLvlMap )      
        aggregate.binaryMap = np.zeros_like( greyLvlMap )        
        aggregate.binaryMap[ np.where( greyLvlMap > aggregate.globalOtsuThreshold ) ] = 1
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e1 = self.calcVoidRatio(aggregate.binaryMap)
        print( 'Void ratio after Otsu binarization= %f' % e1 )
        
        
        # Filling Holes
        aggregate.binaryMap = fillHoles(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        e2 = self.calcVoidRatio(aggregate.binaryMap)        
        print('Void ratio after filling holes = %f' % e2)
        
        
        # Removing specks
        aggregate.binaryMap = removeSpecks(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )        
        e3 = self.calcVoidRatio(aggregate.binaryMap)        
        aggregate.ctVoidRatio = e3       
        print('Void ratio after removing specks = %f' % e3)
        
        # Saving files
        print( "Completed binarization using Otsu threshold of %d" % aggregate.globalOtsuThreshold)                
        otsuBinaryFileName = aggregate.dataOutputDirectory + aggregate.fileName+'-otsuBinary.tiff'
        tiffy.imsave( otsuBinaryFileName, aggregate.binaryMap )        
        otsuTextFileName = aggregate.dataOutputDirectory + aggregate.fileName + '-OtsuThreshold.txt'
        f = open( otsuTextFileName, "w+" )        
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalOtsuThreshold )         

        currentThreshold = aggregate.globalOtsuThreshold
        currentVoidRatio = aggregate.ctVoidRatio
        targetVoidRatio = aggregate.knownVoidRatio
        greyLvlMap = aggregate.filteredGreyLevelMap
            
        tolerance = 0.0001
        deltaVoidRatio = targetVoidRatio - currentVoidRatio
            
        absDeltaThreshold = 0
        maxAbsDeltaThreshold = 500
            
        iterations = 1
        maxIterations = 50
            
        userThresholdStepsTextFileName = aggregate.dataOutputDirectory +  aggregate.fileName + '-userThresholdDensitySteps.txt'
        f = open( userThresholdStepsTextFileName, "w+")
        f.write( "Steps of density based user threshold\n\n" )
        f.close()

        while( abs( deltaVoidRatio ) > tolerance and iterations <= maxIterations ):
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
                
            currentBinaryMap = self.fillholes( currentBinaryMap )
            currentBinaryMap = self.removespecks( currentBinaryMap )
                
            currentVoidRatio = self.calcVoidRatio( currentBinaryMap )
            deltaVoidRatio = targetVoidRatio - currentVoidRatio
                
            print( '\nIteration %d:' % iterations )
            print( 'Current threshold: %d' % round( currentThreshold ) )
            print( 'Otsu Threshold: %d' % round( aggregate.globalOtsuThreshold ) )
            print( 'Current void ratio: %0.5f' % round( currentVoidRatio,5 ) )
            print( 'Target void ratio: %0.5f' % round( targetVoidRatio,5 ) )
            print( 'Target - Current void ratio: %0.5f' % round( deltaVoidRatio, 5 ) )

            f = open( userThresholdStepsTextFileName, 'a' )
            f.write( "\nIteration %d:\n" % iterations)
            f.write( "Current threshold: %d\n" % round( currentThreshold ) )
            f.write( "Otsu Threshold: %d\n" %  round( aggregate.globalOtsuThreshold ) )
            f.write( "Current void ratio:  %0.5f\n" % round( currentVoidRatio,5 ) )
            f.write( "Target void ratio: %0.5f\n" % round( targetVoidRatio,5 ) )
            f.write( "Target - Current void ratio: %0.5f\n" % round( deltaVoidRatio, 5 ) )
            f.close()

            iterations = iterations + 1
            
        aggregate.ctVoidRatio = currentVoidRatio
        aggregate.binaryMapUser = currentBinaryMap
        aggregate.globalUserThreshold = currentThreshold
            
        userBinaryFileName = aggregate.dataOutputDirectory + aggregate.fileName+'-userBinaryDB.tiff'              
        tiffy.imsave( userBinaryFileName, aggregate.binaryMapUser ) 
            
        userThresholdTextFileName = aggregate.dataOutputDirectory +  aggregate.fileName + '-userThresholdDensityBased.txt'
        f = open( userThresholdTextFileName, "w+" )               
        f.write( "\nGlobal User threshold (density based) = %f\n" % aggregate.globalUserThreshold )
        f.close()    
            
        
        
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
        
    
    def fillholes( self, oldBinaryMap ):
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

                   
    def removespecks( self, oldBinaryMap ):
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

        
    def obtainEuclidDistanceMap( self, aggregate ):      
        print('\nFinding Euclidian distance map (EDM)')
        aggregate.euclidDistanceMap = edt( aggregate.binaryMap )             
        print( "EDM Created" )        
        edmImageFileName = aggregate.dataOutputDirectory + aggregate.fileName + '-EDM.tiff'
        tiffy.imsave( edmImageFileName , aggregate.euclidDistanceMap )


    def obtainLocalMaximaMarkers( self, aggregate ):       
        print( '\nFinding local maxima in EDM' )       
        aggregate.markers = localMaxima(aggregate.euclidDistanceMap, allow_borders = False).astype(int)
        
        count=0       
        for i in range( 0, aggregate.markers.shape[ 0 ] ):           
            for j in range( 0, aggregate.markers.shape[ 1 ] ):               
                for k in range( 0, aggregate.markers.shape[ 2 ] ):                   
                    if aggregate.markers[ i ][ j ][ k ] == 1:                       
                        aggregate.markers[ i ][ j ][ k ] = count + 1                       
                        count = count + 1
        makerFileLocationAndName = aggregate.dataOutputDirectory + aggregate.fileName + '-markers.tiff'                
        print( 'Local maxima file created' )
        tiffy.imsave('Markers.tiff',aggregate.markers )


    def obtainLocalHMaximaMarkers( self, aggregate, h ):    
        print( '\nFinding local maxima in EDM' )    
        aggregate.markers = hmax( aggregate.euclidDistanceMap, h ).astype(int)    
        
        count=0    
        for frame in range( 0, aggregate.markers.shape[ 0 ] ):            
            for row in range( 0, aggregate.markers.shape[ 1 ] ):                
                for col in range( 0, aggregate.markers.shape[ 2 ] ):                    
                    if aggregate.markers[ frame ][ row ][ col ] == 1:                        
                        aggregate.markers[ frame ][ row ][ col ] = count + 1                        
                        count = count + 1
        
        makerFileLocationAndName = aggregate.dataOutputDirectory + aggregate.fileName + '-hMarkers.tiff'                
        print( 'Local maxima file created' )
        tiffy.imsave(makerFileLocationAndName, aggregate.markers)


    def obtainTopoWatershed( self, aggregate ):
        print( '\nSegmenting particles by topological watershed' )        
        edm_invert = -aggregate.euclidDistanceMap                  
        binMask = aggregate.binaryMap.astype( bool )              
        particleMarkers = aggregate.markers                    
        aggregate.labelledMap = wsd( edm_invert, markers = particleMarkers, mask = binMask) 
        aggregate.numberOfParticles = aggregate.labelledMap.max()

        for i in range( 1, aggregate.numberOfParticles + 1 ):            
            numberOfParticleVoxel = np.where( aggregate.labelledMap == i)[ 0 ].shape[ 0 ]              
            locData = np.zeros( ( numberOfParticleVoxel, 3 ) )                                        

            for j in range(0, 3):
                locData[ :, j ] = np.where( aggregate.labelledMap == i )[ j ]                            

            p = Particle.Particle( i, numberOfParticleVoxel, locData )                               
            aggregate.particleList.append( p )                                      
        
        print( "Done! List of particles created" )        
        watershedSegmentationFileName = aggregate.dataOutputDirectory + aggregate.fileName + '-watershedSegmentation.tiff'
        tiffy.imsave( watershedSegmentationFileName, aggregate.labelledMap )
        print( 'Watershed segmentation complete' )

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

