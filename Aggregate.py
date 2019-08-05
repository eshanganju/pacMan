# -*- coding: utf-8 -*-

"""
Outline: 
    

Naming conventions:
    Class names: concatenated words each starting with upper case
    Objects, ivars, methods: concatenated words, first word all lower case, subsequent words starting with upper case

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
        self.globalThreshold = 0
        self.binaryMap = np.zeros_like(self.greyLevelMap)
        self.euclidDistanceMap = np.zeros_like(self.binaryMap)
        self.labelledMap = np.zeros_like(self.euclidDistanceMap)
        print("Aggregate activated")
        # TODO: Create and initialize variables as needed in methods