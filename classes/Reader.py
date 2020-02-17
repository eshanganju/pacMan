'''
Reader class
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
        print('\n-----------------')
        print('Reader activated')
        print('-----------------')


    def readTiffStack( self, folderAndFileLocation, cntrZ, cntrY, cntrX, lngt, calib ):       
        print('\nReading tiff stakc from: ' + folderAndFileLocation) 
        upperSlice = cntrZ + round( ( lngt / 2 ) / calib )
        lowerSlice = cntrZ - round( ( lngt / 2 ) / calib )
        upperRow = cntrY + round( ( lngt / 2 ) / calib )
        lowerRow = cntrY - round( ( lngt / 2 ) / calib )
        upperCol = cntrX + round( ( lngt / 2 ) / calib )
        lowerCol = cntrX - round( ( lngt / 2 ) / calib )

        dataBig = tiffy.imread( folderAndFileLocation )
        data = dataBig[ lowerSlice:upperSlice, lowerRow : upperRow, lowerCol : upperCol ] 
        return data 


    def readTiffFileSequence( self, folderLocation, cntrZ, cntrY, cntrX, lngt, calib):
        print( 'Reading tiff files from: ' + folderLocation )

        upperSlice = cntrZ + round( ( lngt / 2 ) / calib )
        lowerSlice = cntrZ - round( ( lngt / 2 ) / calib )
        upperRow = cntrY + round( ( lngt / 2 ) / calib )
        lowerRow = cntrY - round( ( lngt / 2 ) / calib )
        upperCol = cntrX + round( ( lngt / 2 ) / calib )
        lowerCol = cntrX - round( ( lngt / 2 ) / calib )
                
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
