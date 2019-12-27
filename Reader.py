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
from skimage.util import img_as_float
import glob
import matplotlib.pyplot as plt
import gc



# %% Class

class Reader:
    
    '''
    NOTES:
        May need to add an invert fuction as the cone scans are inverted
        
    '''

    def __init__(self):
        print("Reader activated")

    def readTiffStack(self, fileName):
        data = tiffy.imread(fileName)
        return data 

    def readTiffSequence(self, folderLocation, fileRangeMin, fileRangeMax):
        gc.collect()
        print("Reading files from: " + folderLocation)

        searchString = folderLocation+'/*tiff'
        fileList = glob.glob(searchString)
        fileList.sort()
        
        numTiffFiles = (fileRangeMax-fileRangeMin+1)
        tempFile = tiffy.imread(str(fileList[0]))
        rows = tempFile.shape[0]
        columns = tempFile.shape[1]
        del tempFile
        
        gliMap = np.empty((numTiffFiles,rows,columns))
        print('Number of files to read: '+ str(numTiffFiles))    
        
        # Reading individual files - file number offset
        for i in range(fileRangeMin-1,fileRangeMax):
            
            gliMap[ i - fileRangeMin + 1 ] = tiffy.imread ( fileList [ i ] )    
            print('Read ' + str( i - fileRangeMin + 2 ) + '/' + str( fileRangeMax - fileRangeMin + 1 ) + ' files...')
        
        print('\nFinished reading files...')
        self.plotGLI(gliMap)      
        return gliMap
        
    def plotGLI(self,arrayGLI):
        fig, (ax1, ax2, ax3) = plt.subplots(1,3) 
        ax1.imshow(arrayGLI[1], cmap='gray')
        ax2.imshow(arrayGLI[arrayGLI.shape[0]//2], cmap='gray')
        ax3.imshow(arrayGLI[arrayGLI.shape[0]-2], cmap='gray')          
    
    def cropGLI(self, arrayGLI,y1,y2,x1,x2): 
        cropGliMap = arrayGLI[:,y1:y2,x1:x2]
        return cropGliMap
        
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
            


 
