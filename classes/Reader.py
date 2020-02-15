'''

'''


# General
from skimage.util import img_as_float
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt
import numpy as np
import glob
import gc

class Reader:

    def __init__( self ):
        print("Reader activated")


    def readTiffStack( self, fileName ):       
        data = tiffy.imread( fileName )
        return data 


    def readTiffSequence( self, folderLocation, centerSlice, centerRow, centerCol, edgeLength, calib):
        
        plot = LemmeC.LemmeC()
        
        upperSlice = centerSlice + round( ( edgeLength / 2 ) / calib )
        lowerSlice = centerSlice - round( ( edgeLength / 2 ) / calib )
        upperRow = centerRow + round( ( edgeLength / 2 ) / calib )
        lowerRow = centerRow - round( ( edgeLength / 2 ) / calib )
        upperCol = centerCol + round( ( edgeLength / 2 ) / calib )
        lowerCol = centerCol - round( ( edgeLength / 2 ) / calib )
                
        print( 'Reading files from: ' + folderLocation )      
        searchString = folderLocation + '\*tiff'       
        fileList = glob.glob( searchString )        
        fileList.sort()
        
        numTiffFilesInFolder = len(fileList)
        numTiffFilesToRead = ( upperSlice - lowerSlice)
        
        print('Number of files in the folder %d' % numTiffFilesInFolder)
        print('Number of files to read %d' % numTiffFilesToRead)
        
        tempFile = tiffy.imread( str( fileList[ 0 ] ) ) 
        rows = tempFile.shape[ 0 ]
        columns = tempFile.shape[ 1 ]
        del tempFile
        
        gliMap = np.empty( ( numTiffFilesToRead, rows, columns ) ) 
        
        # Reading individual files - file number offset
        for i in range( lowerSlice , upperSlice ):
            gliMap[ i - lowerSlice ] = tiffy.imread ( fileList [ i ] )    
            print('Read ' + str( i - lowerSlice + 1 ) + '/' + str( numTiffFilesToRead ) + ' files...')
        
        print( '\nFinished reading files...' ) 
        croppedGLIMap = gliMap[ :, lowerRow : upperRow, lowerCol : upperCol ] 
        return croppedGLIMap
                
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
            


 
