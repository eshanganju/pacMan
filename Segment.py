# -*- coding: utf-8 -*-

"""
Outline:
This a part of the PAC code; this module is for segmentation of CT data.

---
Features:

---
References: 
    [1] Scikit-image thresholding: https://scikit-image.org/docs/0.13.x/auto_examples/xx_applications/plot_thresholding.html
    [2] Otsu, Nobuyuki. 1979. “A Threshold Selection Method from Gray-Level Histograms.” 
        IEEE Transactions on Systems, Man, and Cybernetics 9 (1): 62–66. 
        https://doi.org/10.1109/TSMC.1979.4310076.

"""

# %% Importing libraries

# Python libraries
import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import local_maxima as localMaxima
from scipy.ndimage.morphology import distance_transform_edt as edt
from skimage.morphology import watershed as wsd

import Particle

#%% Initialize and methods

class Segment:

    # Initialize
    def __init__(self):
        print("Divide-and-conquer: Segmentation tool activated")


    # Binarization - Otsu
    def binarizeOtsu(self,aggregate):
        greyLvlMap=aggregate.greyLevelMap
        aggregate.globalThreshold = threshold_otsu(greyLvlMap)
        aggregate.binaryMap = np.zeros_like(greyLvlMap)
        aggregate.binaryMap[np.where(greyLvlMap >= aggregate.globalThreshold)] = 1
        print("Completed binarization using OTSU")



    # Euclidean Distance map
    def euclidDistanceMap(self,aggregate):
        aggregate.euclidDistanceMap = edt(aggregate.binaryMap)
        print("EDM Created")


    # Location of markers
    def localMaximaMarkers(self,aggregate):
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



    # Topological watershed
    def topoWatershed(self, aggregate):
        
        
        print("Segmenting particles by topological watershed")
        edm_invert=-aggregate.euclidDistanceMap
        binMask = aggregate.binaryMap.astype(bool)
        particleMarkers = aggregate.markers
        aggregate.labelledMap = wsd(edm_invert, markers=particleMarkers, mask=binMask)

        # Creating list of Particles within agregate
        aggregate.particleList=[]
        aggregate.numberOfParticles = aggregate.labelledMap.max()
        aggregate.particleList.append(None)                                         # Done so that list starts from particle index

        for i in range(1,aggregate.numberOfParticles+1):
            numberOfParticleVoxel = np.where(aggregate.labelledMap == i)[0].shape[0]    # Number of voxels also is the volume of the particle
            locData = np.zeros((numberOfParticleVoxel,3))                               # Location data because it is easier to store just the particle loc
            for j in range(0, 3):
                locData[:,j]=np.where(aggregate.labelledMap == i)[j]                    # Extract the locations of the label
            p=Particle.Particle(i,numberOfParticleVoxel,locData)                        # Create particle object
            aggregate.particleList.append(p)                                            # Add particle object to list of particle store in aggregate
        print("Done! List of particles created")



    # Random walker
    def randomWalkerWatershed(self):
        print("Segmentation using random walker watershed")