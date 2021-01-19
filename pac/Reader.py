'''
This library assists in reading and cropping the CT data files
'''

# General imports
from skimage.util import img_as_float
import tifffile as tiffy
import matplotlib.pyplot as plt
import numpy as np
import glob
import gc

# Set this to false if you want to suppress the messages.
VERBOSE = True

def readTiffStack(folderAndFileLocation,cntrZ,cntrY,cntrX,lngt,calib):
    """
    Description:
        Reads a tiff stack created in fiji and passess a subregion of that slice
    Parameters:
        folderAndFileLocation ('str'):
        cntrZ ('int'): Z location of center slice of subregion
        cntrY ('int'): Y location of center slice of subregion
        cntrX ('int'): X location of center slice of subregion
        lngt ('float'): length of subregion in mm
        calib ('float'): calibration factor converting pixel to mm
    Output:
        numpy array of the file.
    """
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

def readDataFromCsv(fileLocation,skipHeader=0,maxRows=0,delim=',',fmt='float',dataForm='number'):
    """
    Description:
        Reads data from csv and converts them to appropriate formats
    Parameters:
        fileLocation (str):
        delim (str):
        skipHeader (int):
        skipFooter (int):
        format (str): int, [float].
        dataForm (str): array or [number]
    Return:
        Data in correct format.
    """
    data = np.genfromtxt(fileLocation,delimiter=delim,skip_header=skipHeader,max_rows=maxRows)

    if dataForm == 'array':
        if fmt == 'float':
            return data.astype('float')
        elif fmt == 'int':
            return data.astype('int')
    else:
        if fmt == 'float':
            return float(data)
        elif fmt == 'int':
            return int(data)

def extractSubregionFromTiffSequence(folderDir,centerZ,topLeftY,topLeftX,lngt,calib,invImg=False,saveImg=False,outputDir='',sampleName='XXX'):
    """
    Description:
        This function reads data from sequence of tiff images from XCT
        The entire image is read and a cubic subregion specified by the user is returned

    Parameters:
        folderDir (str) :   location of the image sequences
        centerZ (int) : slice number of center of the subregion
        topLeftY (int) : Y (vertical) location in pixel units of the lop left corner of the subregion (in correct orientation)
        topLeftX (int) : X (horizontal) location in pixel units of the lop left corner of the subregion (in correct orientation)
        lngt (float) : length in mm of the subregion volume
        calib (float) : calibration factor (mm/pixel)
        invImg (bool) : Keep true if the image needs to be inverted before extraction
        saveImg (bool): Save image in outputDir?
        outputDir (str): path to storage of the filtered image
        sampleName (str): name of the sample used for the analysis, default to sampleX

    Output:
        subregion (np array) : 3D numpy array of extracted subregion.
    """
    # Calculating pixel limits
    upperSlice = centerZ + int( ( lngt / 2 ) / calib)
    lowerSlice = centerZ - int( ( lngt / 2 ) / calib )

    topRow = int(topLeftY)
    bottomRow = int(topLeftY) + int(lngt/calib)

    leftCol = int(topLeftX)
    rightCol = int(topLeftX) + int(lngt/calib)

    print( '\n\nReading files from: ' + folderDir )
    searchString = folderDir + '/*tif'
    searchString2 = folderDir + '/*tiff'

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

    print( 'Number of files in the folder %d' % numTiffFilesInFolder )
    print( 'Number of files to read %d' % numTiffFilesToRead )
    print( '\t' + str( round( lowerSlice ) ) + '-' + str( round(upperSlice) ) )

    tempFile = tiffy.imread( str( fileList[ 0 ] ) )
    rows = tempFile.shape[ 0 ]
    columns = tempFile.shape[ 1 ]
    del tempFile

    gliMap = np.empty( ( numTiffFilesToRead, rows, columns ) )

    # Reading individual files - file number offset
    for i in range( lowerSlice , upperSlice ):
        gliMap[ i - lowerSlice ] = tiffy.imread ( fileList [ i ] )
        print('Read ' + str( i - lowerSlice + 1 ) + '/' + str( numTiffFilesToRead ) + ' files...')

    # Inverting after reading - is this intelligent?
    if invImg == True: gliMap = invertImage(gliMap)
    print( '\nFinished reading files...' )

    # Cropping image
    croppedGLIMap = gliMap[ :, topRow : bottomRow, leftCol : rightCol ]
    if saveImg == True: 
        if VERBOSE: print('\nSaving gli subregion map...')
        tiffy.imsave( outputDir + sampleName + '-gliMap.tif',croppedGLIMap.astype('uint16') )
    return croppedGLIMap

def readTiffFileSequence( folderLocation, cntrZ, cntrY, cntrX, lngt, calib, invImg=False):
    """
    Description:
        This function reads data from sequence of tiff images from XCT
        The entire image is not read, instead a cubic subregion specified by the user is extracted from the scan

    Parameters:
        folderLocation (str) : location of the image sequences
        cntrZ (int) : slice number of center of the subregion
        cntrY (int) : Y (vertical) location in pixel units of the center of the subregion (in correct orientation)
        cntrX (int) : X (horizontal) location in pixel units of the center of the subregion (in correct orientation)
        lngt (float) : length in mm of the subregion volume
        calib (float) : calibration factor (mm/pixel)
        invImg (bool) : Keep true if the image needs to be inverted before extraction

    Output:
        subregion (np array) : 3D numpy array of extracted subregion.
    """
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

    if invImg == True: gliMap = invertImage(gliMap)
    print( '\nFinished reading files...' )
    croppedGLIMap = gliMap[ :, lowerRow : upperRow, lowerCol : upperCol ]
    return croppedGLIMap

def invertImage(gliMapToInvert ):
    """
    Description:

    Parameters:

    Output:
    """
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
