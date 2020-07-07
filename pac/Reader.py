'''
This library assists in reading and cropping the CT data files
'''

# General imports
from skimage.util import img_as_float
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt
import numpy as np
import glob
import gc

def readTiffStack( folderAndFileLocation, cntrZ, cntrY, cntrX, lngt, calib ):
    print('\nReading tiff stack from: ' + folderAndFileLocation)
    upperSlice = cntrZ + round( ( lngt / 2 ) / calib )
    lowerSlice = cntrZ - round( ( lngt / 2 ) / calib )
    upperRow = cntrY + round( ( lngt / 2 ) / calib )
    lowerRow = cntrY - round( ( lngt / 2 ) / calib )
    upperCol = cntrX + round( ( lngt / 2 ) / calib )
    lowerCol = cntrX - round( ( lngt / 2 ) / calib )

    dataBig = tiffy.imread( folderAndFileLocation )
    data = dataBig[ lowerSlice:upperSlice, lowerRow : upperRow, lowerCol : upperCol ]
    return data

def readTiffFileSequence( folderLocation, cntrZ, cntrY, cntrX, lngt, calib, invertImageData=True):
    print( 'Reading tiff files from: ' + folderLocation )

    upperSlice = cntrZ + round( ( lngt / 2 ) / calib )
    lowerSlice = cntrZ - round( ( lngt / 2 ) / calib )
    upperRow = cntrY + round( ( lngt / 2 ) / calib )
    lowerRow = cntrY - round( ( lngt / 2 ) / calib )
    upperCol = cntrX + round( ( lngt / 2 ) / calib )
    lowerCol = cntrX - round( ( lngt / 2 ) / calib )

    print( 'Reading files from: ' + folderLocation )
    searchString = folderLocation + '/*tif'
    searchString2 = folderLocation + '/*tiff'

    fileList1 = glob.glob(searchString)
    fileList2 = glob.glob(searchString2)

    if len(fileList1)>len(fileList2):
        fileList = fileList1
        fileList.sort()
    else:
        fileList = fileList2
        fileList.sort()

    numTiffFilesInFolder = len(fileList)
    numTiffFilesToRead = ( upperSlice - lowerSlice)

    print('Number of files in the folder %d' % numTiffFilesInFolder)
    print('Number of files to read %d' % numTiffFilesToRead)
    print(str(lowerSlice) + '-' + str(upperSlice))

    tempFile = tiffy.imread( str( fileList[ 0 ] ) )
    rows = tempFile.shape[ 0 ]
    columns = tempFile.shape[ 1 ]
    del tempFile

    gliMap = np.empty( ( numTiffFilesToRead, rows, columns ) )

    # Reading individual files - file number offset
    for i in range( lowerSlice , upperSlice ):
        gliMap[ i - lowerSlice ] = tiffy.imread ( fileList [ i ] )
        print('Read ' + str( i - lowerSlice + 1 ) + '/' + str( numTiffFilesToRead ) + ' files...')

    if invertImageData == True: gliMap = invertImage(gliMap)
    print( '\nFinished reading files...' )
    croppedGLIMap = gliMap[ :, lowerRow : upperRow, lowerCol : upperCol ]
    return croppedGLIMap

def invertImage( gliMapToInvert ):
    invertedGliMap = gliMapToInvert

    invertZ = input('\nInvert Z direction(y/[n]):')
    if invertZ.lower() == 'y': invertZ = True
    else: invertZ = False

    invertY = input('Invert Y direction(y/[n]):')
    if invertY.lower() == 'y': invertY = True
    else: invertY = False

    invertX = input('Invert X direction(y/[n]):')
    if invertX.lower() == 'y': invertX = True
    else: invertX = False

    if invertZ == True:
        invertedGliMap = np.flip( invertedGliMap , 0 )

    if invertY == True:
        invertedGliMap = np.flip( invertedGliMap , 1 )

    if invertX == True:
        invertedGliMap = np.flip( invertedGliMap , 2 )

    return invertedGliMap
