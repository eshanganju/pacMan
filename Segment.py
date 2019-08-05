# -*- coding: utf-8 -*-

"""

References: 
    [1] Scikit-image thresholding: https://scikit-image.org/docs/0.13.x/auto_examples/xx_applications/plot_thresholding.html

"""

# %% Importing libraries

# Python libraries
import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import local_maxima as localMaxima
from scipy.ndimage.morphology import distance_transform_edt as edt
from skimage.morphology import watershed as wsd

#%% Initialize and methods

class Segment:

    # Initialize
    def __init__(self):
        print("Segmenter activated")



    # Binarization - Otsu
    def binarizeOtsu(self,greyLvlMap):
        '''
        Binarization converts grey level intensity map into a binary map
        Scikit-image libraries have some good options [1]
        We stick to Otsu's method [2] for now
        '''
        print("Starting binarization using OTSU")
        globalThresholdIntensity = threshold_otsu(greyLvlMap)
        binaryMap = np.zeros_like(greyLvlMap)
        binaryMap[np.where(greyLvlMap >= globalThresholdIntensity)] = 1
        return globalThresholdIntensity, binaryMap



    # Euclidean Distance map
    def euclidDistanceMap(self,binaryMap):
        print("Creating EDM")
        edm = edt(binaryMap)
        return edm



    # Topological watershed
    def topoWatershed(self, binaryMap, distanceMap):
        print("Segmenting particles by topological watershed")
        



    # Random walker
    def randomWalkerWatershed(self):
        print("Segmentation using random walker watershed")