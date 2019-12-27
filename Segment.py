# -*- coding: utf-8 -*-

"""
Outline:
This a part of the PAC code; this module is for segmentation of CT data.

---
Features:

---
References: 


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

import Particle

import skimage.external.tifffile as tiffy

import matplotlib.pyplot as plt


# %% Class

class Segment:

    # Initialize
    def __init__( self ):
        print( "Segmenter activated" )


    # Binarization - Otsu
    def binarizeOtsu( self, aggregate ):
        
        greyLvlMap = aggregate.greyLevelMap
        
        aggregate.globalOtsuThreshold = threshold_otsu( greyLvlMap )
        
        aggregate.binaryMap[ np.where( greyLvlMap > aggregate.globalOtsuThreshold ) ] = 1
        
        aggregate.binaryMap = aggregate.binaryMap.astype(int)
        
        print( "\n\nCompleted binarization using OTSU" )
        
        tiffy.imsave( 'Bin.tiff', aggregate.binaryMap )
        
        f = open( "OtsuThreshold.txt", "w+" )
        f.write( "\nGlobal Otsu threshold = %f\n" % aggregate.globalOtsuThreshold )
        f.close()

        
    # Fill holes
    def fillholes( self, aggregate ):
        
        oldBinaryMap = aggregate.binaryMap
        
        newBinaryMap = fillHoles( aggregate.binaryMap )
            
        fig1, ( ( ax11, ax12, ax13 ), ( ax21, ax22, ax23 ) ) = plt.subplots( 2, 3 ) 
            
        ax11.imshow( newBinaryMap[ newBinaryMap.shape[ 0 ] // 4 ], cmap = 'gray' )
            
        ax12.imshow( newBinaryMap[ newBinaryMap.shape[ 0 ] // 2 ], cmap = 'gray' )
            
        ax13.imshow( newBinaryMap[ ( newBinaryMap.shape[ 0 ] // 4 ) * 3 ], cmap='gray')  
            
        ax21.imshow( oldBinaryMap[ oldBinaryMap.shape[ 0 ] // 4 ], cmap = 'gray' )
            
        ax22.imshow( oldBinaryMap[ oldBinaryMap.shape[ 0 ] // 2 ], cmap = 'gray' )
            
        ax23.imshow( oldBinaryMap[ ( oldBinaryMap.shape[ 0 ] // 4 ) * 3 ], cmap='gray')            
                   
        aggregate.binaryMap = newBinaryMap
        
    
    def removespecks( self, aggregate ):
        
        oldBinaryMap = aggregate.binaryMap
        
        newBinaryMap = removeSpecks( oldBinaryMap )
        
        fig1, ( ( ax11, ax12, ax13 ), ( ax21, ax22, ax23 ) ) = plt.subplots( 2, 3 ) 
            
        ax11.imshow( newBinaryMap[ newBinaryMap.shape[ 0 ] // 4 ], cmap = 'gray' )
            
        ax12.imshow( newBinaryMap[ newBinaryMap.shape[ 0 ] // 2 ], cmap = 'gray' )
            
        ax13.imshow( newBinaryMap[ ( newBinaryMap.shape[ 0 ] // 4 ) * 3 ], cmap='gray')  
            
        ax21.imshow( oldBinaryMap[ oldBinaryMap.shape[ 0 ] // 4 ], cmap = 'gray' )
            
        ax22.imshow( oldBinaryMap[ oldBinaryMap.shape[ 0 ] // 2 ], cmap = 'gray' )
            
        ax23.imshow( oldBinaryMap[ ( oldBinaryMap.shape[ 0 ] // 4 ) * 3 ], cmap='gray') 
        
        aggregate.binaryMap = newBinaryMap
        


    # Euclidean Distance map
    def euclidDistanceMap( self, aggregate ):
        
        aggregate.euclidDistanceMap = edt( aggregate.binaryMap )
        
        edm = aggregate.euclidDistanceMap
        
        fig1, ( ax11, ax12, ax13 ) = plt.subplots( 1, 3 ) 
            
        ax11.imshow( edm[ edm.shape[ 0 ] // 4 ], cmap = 'gray' )
            
        ax12.imshow( edm[ edm.shape[ 0 ] // 2 ], cmap = 'gray' )
            
        ax13.imshow( edm[ ( edm.shape[ 0 ] // 4 ) * 3 ], cmap='gray')  
        
        print( "EDM Created" )
        
        tiffy.imsave( 'EDM.tiff', aggregate.euclidDistanceMap )


    # Location of markers
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
        
    # Location of markers
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


    # Topological watershed
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
        
        tiffy.imsave( 'watershedSegmentation.tiff', aggregate.labelledMap )

    # Random walker
    def randomWalkerWatershed2Particles( self ):
        print( "Segmentation using random walker watershed for 2 particles" )
