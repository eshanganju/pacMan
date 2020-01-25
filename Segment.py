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
        e = (v-vs)/vs        
        print('void ratio after removing specks = %f' % e)
        
        
        # Saving files
        print( "Completed binarization using Otsu" )                
        tiffy.imsave( 'BinOtsu.tiff', aggregate.binaryMap.astype(int) )        
        f = open( "OtsuThreshold.txt", "w+" )        
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalOtsuThreshold )        
        f.close()
        
        

    # Reset binarization from Otsu to fix diparity with measured void ratio 
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
        e = (v-vs)/vs        
        print('void ratio after removing specks = %f' % e)
        
        
        # Saving files
        print( "Completed binarization using user input" )                
        tiffy.imsave( 'BinUser.tiff', aggregate.binaryMap )        
        f = open( "UserThreshold.txt", "w+" )        
        f.write( "\nGlobal User threshold = %f\n" % aggregate.globalOtsuThreshold )        
        f.close()
        
        
   
    # Fill holes
    def fillholes( self, oldBinaryMap ):
        newBinaryMap = fillHoles( oldBinaryMap )
        return newBinaryMap.astype(int) 
                
    # Rem specks    
    def removespecks( self, oldBinaryMap ):
        newBinaryMap = removeSpecks( oldBinaryMap )   
        return newBinaryMap.astype(int)
        
        
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
        tiffy.imsave( 'watershedSegmentation.tif', aggregate.labelledMap )
   
     
    # Reset particleList
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

   
    # Random walker
    def randomWalkerWatershed2Particles( self ):
        print( "Segmentation using random walker watershed for 2 particles" )
