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
from scipy.ndimage.morphology import distance_transform_edt as edt
from skimage.morphology import watershed as wsd
import matplotlib.pyplot as plt
import Particle
import skimage.external.tifffile as tiffy

class Segment:

    # Initialize
    def __init__(self):
        print("Segmenter activated")


    # Histogram
    def greyLevelHistogram(self, aggregate):

        # Unfiltered histogram
        greyLevelMap = aggregate.greyLevelMap
        greyLevelList = np.zeros(greyLevelMap.size)
        count = 0
        for i in range(0,greyLevelMap.shape[0]):
            for j in range(0,greyLevelMap.shape[1]):
                for k in range(0,greyLevelMap.shape[2]):
                    greyLevelList[count]=greyLevelMap[i][j][k]
                    count=count+1
        greyLevelHist = np.histogram(greyLevelList, bins=1000, range=(0,1))
        aggregate.greyLevelHistogram[:,0] = np.arange(0.0005,1,0.001)
        aggregate.greyLevelHistogram[:,1] = (greyLevelHist[0]/(greyLevelHist[0].sum()))*100

        # Filtered histogram
        filteredGreyLevelMap = aggregate.filteredGreyLevelMap
        filteredGreyLevelList = np.zeros(filteredGreyLevelMap.size)
        count = 0
        for i in range(0,filteredGreyLevelMap.shape[0]):
            for j in range(0,filteredGreyLevelMap.shape[1]):
                for k in range(0,filteredGreyLevelMap.shape[2]):
                    filteredGreyLevelList[count]=filteredGreyLevelMap[i][j][k]
                    count=count+1
        filteredGreyLevelHist = np.histogram(filteredGreyLevelList, bins=1000, range=(0,1))
        aggregate.filteredGreyLevelHistogram[:,0] = np.arange(0.0005,1,0.001)
        aggregate.filteredGreyLevelHistogram[:,1] = (filteredGreyLevelHist[0]/(filteredGreyLevelHist[0].sum()))*100

        # TODO: Get rid of this after visualization is updated
        plt.figure()
        plt.plot(aggregate.greyLevelHistogram[:,0], aggregate.greyLevelHistogram[:,1],'k--',label='Unfiltered image')
        plt.plot(aggregate.filteredGreyLevelHistogram[:,0], aggregate.filteredGreyLevelHistogram[:,1],'k',label='Filtered image')
        plt.legend(loc='upper right')
        plt.xlim(0,1)
        plt.xticks([0,.10,.20,0.25,.30,.40,0.49,.50,.60,.70,0.75,.80,.90,.100])
        plt.ylim(0,5)
        plt.show()


    # Binarization - Otsu
    def binarizeOtsu(self,aggregate):
        greyLvlMap=aggregate.greyLevelMap
        aggregate.globalOtsuThreshold = threshold_otsu(greyLvlMap)
        aggregate.binaryMap[np.where(greyLvlMap > aggregate.globalOtsuThreshold)] = 1
        print("\n\nCompleted binarization using OTSU")
        tiffy.imsave('Bin.tiff',aggregate.binaryMap)
        f = open("OtsuThreshold.txt","w+")
        f.write("\nGlobal Otsu threshold = %f\n" % aggregate.globalOtsuThreshold)
        f.close()



    # Euclidean Distance map
    def euclidDistanceMap(self,aggregate):
        aggregate.euclidDistanceMap = edt(aggregate.binaryMap)
        print("EDM Created")
        tiffy.imsave('EDM.tiff',aggregate.euclidDistanceMap)

    # Location of markers
    def localMaximaMarkers(self,aggregate):
        #TODO: Check implementation if h-maxima for this steps
        print("Finding local maxima in EDM")
        aggregate.markers = localMaxima(aggregate.euclidDistanceMap).astype(int)
        
        # Give each local maxima a new label
        count=0
        for i in range(0,aggregate.markers.shape[0]):
            for j in range(0,aggregate.markers.shape[1]):
                for k in range(0,aggregate.markers.shape[2]):
                    if aggregate.markers[i][j][k] == 1:
                        aggregate.markers[i][j][k] = count + 1
                        count = count + 1
        tiffy.imsave('Markers.tiff',aggregate.markers)

    # Topological watershed
    def topoWatershed(self, aggregate):
        print("Segmenting particles by topological watershed")
        
        edm_invert=-aggregate.euclidDistanceMap          
        binMask = aggregate.binaryMap.astype(bool)      
        particleMarkers = aggregate.markers            
        aggregate.labelledMap = wsd(edm_invert, markers=particleMarkers, mask=binMask) 

        # Creating list of Particles within agregate
        aggregate.numberOfParticles = aggregate.labelledMap.max()

        for i in range(1,aggregate.numberOfParticles+1):
            numberOfParticleVoxel = np.where(aggregate.labelledMap == i)[0].shape[0]  
            locData = np.zeros((numberOfParticleVoxel,3))                            
            for j in range(0, 3):
                locData[:,j]=np.where(aggregate.labelledMap == i)[j]                
            p=Particle.Particle(i,numberOfParticleVoxel,locData)                   
            aggregate.particleList.append(p)                                      
        print("Done! List of particles created")
        tiffy.imsave('watershedSegmentation.tiff', aggregate.labelledMap)

    # Random walker
    def randomWalkerWatershed2Particles(self):
        print("Segmentation using random walker watershed for 2 particles")
