# -*- coding: utf-8 -*-

"""
Outline: 

---
Features:

---
References [All reference will be cited in the method or class]:
    

"""

# %% Import libraries

# Python libraries
import numpy as np

# %% Aggregate class


class Aggregate:

    def __init__(self, greyLvlMap, pxlSize):
        self.pixelSize = pxlSize
        self.greyLevelMap = greyLvlMap
        self.filteredGreyLevelMap = np.zeros_like(self.greyLevelMap)
        self.binaryMap = np.zeros_like(self.greyLevelMap)
        self.globalOtsuThreshold = 0
        self.euclidDistanceMap = np.zeros_like(self.greyLevelMap)
        self.markers = np.zeros_like(self.greyLevelMap)
        self.labelledMap = np.zeros_like(self.greyLevelMap)
        self.particleList = []         
        self.particleList.append(None)  
        self.numberOfParticles = 0      
        self.particleSizeDataSummary = np.zeros((5000,6)) 
        self.grainSizeDistributionEquivalentSphere = np.zeros((5000,2))     
        self.grainSizeDistributionFeretMax = np.zeros((5000,2))         
        self.grainSizeDistributionFeretMin = np.zeros((5000,2))             
        self.grainSizeDistributionFeretMed = np.zeros((5000,2))             
        self.greyLevelHistogram = np.zeros((1000,2))                    
        self.benchMarkNumberOfParticles = 0
        self.benchMarkCentres = np.zeros((self.benchMarkNumberOfParticles,3))
        self.benchMarkRadii = np.zeros((self.benchMarkNumberOfParticles,1))
        self.benchMarkGrainSizeDistribution = np.zeros((self.benchMarkNumberOfParticles,2))
        self.benchMarkNumberOfContacts = 0
        self.benchMarkContactNormal = np.zeros((1,5))

        print("Aggregate activated")


