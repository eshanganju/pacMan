# -*- coding: utf-8 -*-

"""
Purpose: 
    This class will read the data files from the local storage
    Convert the GLI from what is on the data to 0-1 range
    Crop image to get appropriate section of the data
    
"""

"""
State: 
    Normalizing large arrays - how?
"""

# General
import numpy as np
import skimage.external.tifffile as tiffy
import glob
import matplotlib.pyplot as plt


# %% Class

class Reader:

    def __init__(self):
        print("Reader activated")

    def readTiffStack(self, fileName):
        data = tiffy.imread(fileName)
        return data 

    def readTiffSequence(self, folderLocation, fileExtension):
        print("Reading files from: " + folderLocation)

        # Reading files
        searchString = folderLocation+'/*'+fileExtension
        numTiffFiles = len(glob.glob(searchString))
        print('Number of files to read: '+ str(numTiffFiles))    
        gliMap = tiffy.imread(searchString)  
        print('Finished reading files...')
        
        # Normalization
        print('/nStarting normalization...')
        minVal=0
        maxVal=self.maxGliPossible(gliMap)
        normGliMap=self.normalizeArray(gliMap,minVal,maxVal)
        
        print('\nFiles read and normalized...')
        return normGliMap

    def normalizeArray(self, arrayGLI, minVal, maxVal):
        print("Normalizing array...\n")
        
        minVal = float(minVal)
        maxVal = float(maxVal)
        
        normArrayGli = np.zeros_like(arrayGLI)
        
        for i in range(0,arrayGLI.shape[0]):
            normArrayGli[i]=(float(arrayGLI[i])-minVal)/(maxVal-minVal)
        return normArrayGli

    def maxGliPossible(self, arrayGli):
        dataType = str(arrayGli.dtype)   
         
        if dataType=='uint8':
            maxValOption = (2**8-1)

        elif dataType=='uint16':
            maxValOption = (2**16-1)
        
        elif dataType=='uint32':
            maxValOption = (2**32-1)
        
        elif dataType=='uint64':
            maxValOption = (2**64-1)
        
        return maxValOption
            


 
