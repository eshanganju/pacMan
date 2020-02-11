# -*- coding: utf-8 -*-

"""
Outline:
This a part of the PAC code; this module is for segmentation of CT data.
---
Features:
    Binarize
    Fix errors - fill holes and remove specks in binary image
    Find Euclidian distance map (EDM)
    Locate maxima in EDM
        local maxima
        local-h-maxima
    Topoglogical watershed segmentation
    Fix particle list
    Random walker segmentation 2 particles
---
References: 
    
---
TODO: 
    Binarization - how to ensure particle size correction?


"""

# Importing libraries
import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import local_maxima as localMaxima
from skimage.morphology import h_maxima as hmax
from scipy.ndimage.morphology import distance_transform_edt as edt
from scipy.ndimage.morphology import binary_fill_holes as fillHoles
from scipy.ndimage.morphology import binary_opening as removeSpecks
from skimage.morphology import watershed as wsd
import classes.Particle as particle
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt


# %% Class

class Segment:

    def __init__( self ):
        print( "Segmenter activated" )


    def binarizeOtsu( self, aggregate ):
        
        greyLvlMap = aggregate.filteredGreyLevelMap
        aggregate.globalOtsuThreshold = threshold_otsu( greyLvlMap )      
        aggregate.binaryMap = np.zeros_like( greyLvlMap )        
        aggregate.binaryMap[ np.where( greyLvlMap > aggregate.globalOtsuThreshold ) ] = 1
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        vs = aggregate.binaryMap.sum()
        v = ( aggregate.binaryMap.shape[ 0 ]) * aggregate.binaryMap.shape[ 1 ] * aggregate.binaryMap.shape[ 2 ]
        e = ( v - vs ) / vs
        print( 'void ratio after Otsu binarization= %f' % e )
        
        
        # Filling Holes
        aggregate.binaryMap = fillHoles(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        vs = aggregate.binaryMap.sum()      
        v = (aggregate.binaryMap.shape[0])*aggregate.binaryMap.shape[1]*aggregate.binaryMap.shape[2]        
        e = (v-vs)/vs        
        print('void ratio after filling holes = %f' % e)
        
        
        # Removing specks
        aggregate.binaryMap = removeSpecks(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )        
        vs = aggregate.binaryMap.sum()       
        v = (aggregate.binaryMap.shape[0])*aggregate.binaryMap.shape[1]*aggregate.binaryMap.shape[2]        
        aggregate.ctVoidRatio = (v-vs)/vs        
        print('void ratio after removing specks = %f' % e)
        
        # Saving files
        print( "Completed binarization using Otsu threshold of %d" % aggregate.globalOtsuThreshold)                
        otsuBinaryFileName = aggregate.fileName+'-otsuBinary.tiff'
        tiffy.imsave( otsuBinaryFileName, aggregate.binaryMap )        
        otsuTextFileName = aggregate.fileName + '-OtsuThreshold.txt'
        f = open( otsuTextFileName, "w+" )        
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalOtsuThreshold )        
        f.close()

               
    def resetOtsuBinarizationAccordingToUser( self, aggregate, userThreshold):
        
        # Resetting binarization threshold
        greyLvlMap = aggregate.filteredGreyLevelMap       
        aggregate.globalUserThreshold = userThreshold        
        aggregate.binaryMap = np.zeros_like( greyLvlMap )        
        aggregate.binaryMap[ np.where( greyLvlMap > aggregate.globalUserThreshold ) ] = 1
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        vs = aggregate.binaryMap.sum()
        v = aggregate.binaryMap.shape[ 0 ] * aggregate.binaryMap.shape[ 1 ] * aggregate.binaryMap.shape[ 2 ]
        e = ( v - vs ) / vs
        print( 'void ratio after binarization= %f' % e )
                
        # Filling Holes
        aggregate.binaryMap = fillHoles(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        vs = aggregate.binaryMap.sum()      
        v = aggregate.binaryMap.shape[0]*aggregate.binaryMap.shape[1]*aggregate.binaryMap.shape[2]        
        e = (v-vs)/vs        
        print('void ratio after filling holes = %f' % e)
                
        # Removing specks
        aggregate.binaryMap = removeSpecks(aggregate.binaryMap)
        aggregate.binaryMap = aggregate.binaryMap.astype( int )
        vs = aggregate.binaryMap.sum()       
        v = aggregate.binaryMap.shape[0]*aggregate.binaryMap.shape[1]*aggregate.binaryMap.shape[2]        
        aggregate.ctVoidRatio = (v-vs)/vs        
        print('void ratio after removing specks = %f' % e)
                
        # Saving files
        print( "Completed binarization using user input" )  
        userBinaryFileName = aggregate.fileName+'-userBinary.tiff'              
        tiffy.imsave( userBinaryFileName, aggregate.binaryMapUser )        
        
        userThresholdTextFileName = aggregate.fileName + '-userThreshold.txt'
        f = open( userThresholdTextFileName, "w+" )               
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalUserThreshold )        
        f.close()        

    def resetBinarizationAccordingToDensity( self, aggregate ):
        '''
        
        Parameters
        ----------
        aggregate : Object
            Object of the sand particles containing atleast the GLI
        
        targetVoidRatio : float
            Void ratio measured from 

        Returns
        -------
        Binary map according to void ratio measured from samples

        '''
        if aggregate.globalOtsuThreshold != 0:
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
                
                iterations = iterations + 1
            
            aggregate.ctVoidRatio = currentVoidRatio
            aggregate.binaryMapUser = currentBinaryMap
            aggregate.globalUserThreshold = currentThreshold
            
            userBinaryFileName = aggregate.fileName+'-userBinaryDB.tiff'              
            tiffy.imsave( userBinaryFileName, aggregate.binaryMapUser ) 
            
            userThresholdTextFileName = aggregate.fileName + '-userThresholdDensityBased.txt'
            f = open( userThresholdTextFileName, "w+" )               
            f.write( "\nGlobal User threshold (density based) = %f\n" % aggregate.globalUserThreshold )        
            f.close()
            
        else:
            print('Run Otsu first!')
        
        
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
            DESCRIPTION.

        '''
        newBinaryMap = fillHoles( oldBinaryMap )
        return newBinaryMap.astype(int) 

                   
    def removespecks( self, oldBinaryMap ):
        newBinaryMap = removeSpecks( oldBinaryMap )   
        return newBinaryMap.astype(int)

        
    def euclidDistanceMap( self, aggregate ):      
        aggregate.euclidDistanceMap = edt( aggregate.binaryMap )
        edm = aggregate.euclidDistanceMap        
        
        fig1, ( ax11, ax12, ax13 ) = plt.subplots( 1, 3 )             
        ax11.imshow( edm[ edm.shape[ 0 ] // 4 ], cmap = 'gray' )            
        ax12.imshow( edm[ edm.shape[ 0 ] // 2 ], cmap = 'gray' )            
        ax13.imshow( edm[ ( edm.shape[ 0 ] // 4 ) * 3 ], cmap='gray')          
        
        print( "EDM Created" )        
        tiffy.imsave( 'EDM.tiff', aggregate.euclidDistanceMap )


    def localMaximaMarkers( self, aggregate ):       
        print( "Finding local maxima in EDM" )       
        aggregate.markers = localMaxima(aggregate.euclidDistanceMap, allow_borders = False).astype(int)
        
        count=0       
        for i in range( 0, aggregate.markers.shape[ 0 ] ):           
            for j in range( 0, aggregate.markers.shape[ 1 ] ):               
                for k in range( 0, aggregate.markers.shape[ 2 ] ):                   
                    if aggregate.markers[ i ][ j ][ k ] == 1:                       
                        aggregate.markers[ i ][ j ][ k ] = count + 1                       
                        count = count + 1
                        
        tiffy.imsave('Markers.tiff',aggregate.markers )


    def localhMaxima( self, aggregate, h ):    
        
        print( "Finding local maxima in EDM" )    
        aggregate.markers = hmax( aggregate.euclidDistanceMap, h ).astype(int)    
        
        count=0    
        for i in range( 0, aggregate.markers.shape[ 0 ] ):            
            for j in range( 0, aggregate.markers.shape[ 1 ] ):                
                for k in range( 0, aggregate.markers.shape[ 2 ] ):                    
                    if aggregate.markers[ i ][ j ][ k ] == 1:                        
                        aggregate.markers[ i ][ j ][ k ] = count + 1                        
                        count = count + 1
                        
        tiffy.imsave('Markers.tiff',aggregate.markers)


    def topoWatershed( self, aggregate ):
        
        print( "Segmenting particles by topological watershed" )        
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
        tiffy.imsave( 'watershedSegmentation.tif', aggregate.labelledMap )

    def fixErrorsInSegmentation(self, aggregate):
        '''
        Parameters
        ----------
        aggregate : Object
            

        Returns
        -------
        Updates aggregate object to correct oversegmentation and fixes the particle lists

        '''
        
    def resetParticleList(self, aggregate):
        
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


    def randomWalkerWatershed2Particles( self ):
        print( "Segmentation using random walker watershed for 2 particles" )
