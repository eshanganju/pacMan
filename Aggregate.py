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
        self.binaryMap = np.zeros_like(self.greyLevelMap)
        self.globalThreshold = 0
        self.euclidDistanceMap = np.zeros_like(self.greyLevelMap)
        self.markers = np.zeros_like(self.greyLevelMap)
        self.labelledMap = np.zeros_like(self.greyLevelMap)
        self.particleList = []                                              # TODO: Object of objects? List of objects is irritating, or not?
        self.particleList.append(None)                                      # Done to start at index 1, which matches with particle index (0 is void)
        self.numberOfParticles = 0
        self.grainSizeDistributionEquivalentSphere = np.zeros((5000,2))     # GSD sphere (Size is arbitrary - rows will correspond to number of particles)
        self.grainSizeDistributionFeretMax = np.zeros((5000,2))             # GSD FeretMax ([0]Size is arbitrary - rows = numberOfParticles
        self.grainSizeDistributionFeretMin = np.zeros((5000,2))             # GSD FeretMin ([0]Size is arbitrary - rows = numberOfParticles
        self.grainSizeDistributionFeretMed = np.zeros((5000,2))             # GSD FeretMed ([0]Size is arbitrary - rows = numberOfParticles
        self.greyLevelHistogram = np.zeros((1000,2))                        # Histogram with x axis (0-1)
        print("Its like a beach in here: Aggregate activated")


