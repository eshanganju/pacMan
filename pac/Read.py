"""This library assists in reading and cropping the CT data files
"""

# General imports
from skimage.util import img_as_float
import tifffile as tiffy
import matplotlib.pyplot as plt
import numpy as np
import glob
import gc

# Set this to false if you want to suppress the messages.
TESTING = True
VERBOSE = True

def readDataFromCsv(fileLocation,skipHeader=0,maxRows=0,delim=',',fmt='float',dataForm='number'):
    """Reads data from csv and converts them to appropriate formats

    Parameters
    ----------
    fileLocation : str
    delim : str
    skipHeader : int
    skipFooter : int
    fmt : str
        int or [float]
    dataForm : str
        array or [number]

    Return
    ------
    data : number or ndarray in int or float
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

def extractSubregionFromTiffSequence( folderDir, Z, Y, X, lngt, calib, reference='topLeft',
                                      invImg=False, saveImg=False, outputDir='', sampleName='XXX'):
    """This function reads data from sequence of tiff images from XCT
    The entire image is read and a cubic subregion specified by the user is returned

    Parameters
    ----------
    folderDir : str
        location of the image sequences
    reference : str
        topLeft - the ZYX locations are of the top left
        center - the ZYX locations are of the the center
    Z : int
        slice number of center of the subregion
    Y : int
        Y (vertical) location in pixel units of the subregion
         (in correct orientation)
    X : int
        X (horizontal) location in pixel units of the subregion
        (in correct orientation)
    lngt : float
        length in mm of the subregion volume
    calib : float
        calibration factor (mm/pixel)
    invImg : bool
        Keep true if the image needs to be inverted before extraction
    saveImg : bool
        Save image in outputDir?
    outputDir : str
        Path to storage of the filtered image
    sampleName : str
        name of the sample used for the analysis, default to sampleX

    Output
    ------
    subregion : unsigned integer ndarray
        3D numpy array of extracted subregion.
    """
    upperSlice = Z + int( ( lngt / 2 ) / calib)
    lowerSlice = Z - int( ( lngt / 2 ) / calib )

    if reference == 'topLeft':
        topRow = int(Y)
        bottomRow = int(Y) + int(lngt/calib)
        leftCol = int(X)
        rightCol = int(X) + int(lngt/calib)

    elif reference == 'center':
        topRow = int(Y) - int(lngt/calib)//2
        bottomRow = int(Y) + int(lngt/calib)//2
        leftCol = int(X) - int(lngt/calib)//2
        rightCol = int(X) + int(lngt/calib)//2

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

    for i in range( lowerSlice , upperSlice ):
        gliMap[ i - lowerSlice ] = tiffy.imread ( fileList [ i ] )
        print('Read ' + str( i - lowerSlice + 1 ) + '/' + str( numTiffFilesToRead ) + ' files...')

    if invImg == True: gliMap = invertImage(gliMap)
    print( '\nFinished reading files...' )

    croppedGLIMap = gliMap[ :, topRow : bottomRow, leftCol : rightCol ]
    if saveImg == True: 
        if VERBOSE: print('\nSaving gli subregion map...')
        tiffy.imsave( outputDir + sampleName + '-gliMap.tif',croppedGLIMap.astype('uint16') )

    return croppedGLIMap

def extractSubregionFromTiffFile( folderDir, 
									Z, Y, X, lngt, 
									calib, 
									zReference = 'middle', xyreference='topLeft', 
									invImg=False, 
									saveImg=False, 
									outputDir='', 
									sampleName='XXX' ):
    """This function reads data from a single tiff stack
    The entire image is read and a cubic subregion specified by the user is returned

    Parameter
    ----------
    folderDir : str
        location of the image
    zReference : str
        center - the Z location is the center
		low - the Z location is the lower Z value 
		high - the Z location is the higher Z value 
    xyReference : str
        topLeft - the YX locations are of the top left
        center - the YX locations are of the the center
    Z : int
        slice number of the subregion (start is 0)
    Y : int
        Y (vertical) location in pixel units of the subregion
        (in correct orientation)
    X : int
        X (horizontal) location in pixel units of the subregion
        (in correct orientation)
    lngt : float
        length in mm of the subregion volume
    calib : float
        calibration factor (mm/pixel)
    invImg : bool
        Keep true if the image needs to be inverted before extraction
    saveImg : bool
        Save image in outputDir?
    outputDir : str
        Path to storage of the filtered image
    sampleName : str
        name of the sample used for the analysis, default to sampleX

    Output
    ------
    subregion : unsigned integer ndarray
        3D numpy array of extracted subregion.
    """
    if zReference == 'center':
		upperSlice = Z + int( ( lngt / 2 ) / calib )
    	lowerSlice = Z - int( ( lngt / 2 ) / calib )

	elif zReference == 'low':
		upperSlice = Z + int( lngt / calib ) 
		lowerslice = Z 

	elif zReference == 'high':
		upperSlice = Z
		lowerSlice = Z - int( lngt / calib )

    if reference == 'topLeft':
        topRow = int(Y)
        bottomRow = int(Y) + int(lngt/calib)
        leftCol = int(X)
        rightCol = int(X) + int(lngt/calib)

    elif reference == 'center':
        topRow = int(Y) - int(lngt/calib)//2
        bottomRow = int(Y) + int(lngt/calib)//2
        leftCol = int(X) - int(lngt/calib)//2
        rightCol = int(X) + int(lngt/calib)//2

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

    for i in range( lowerSlice , upperSlice ):
        gliMap[ i - lowerSlice ] = tiffy.imread ( fileList [ i ] )
        print('Read ' + str( i - lowerSlice + 1 ) + '/' + str( numTiffFilesToRead ) + ' files...')

    if invImg == True: gliMap = invertImage(gliMap)
    print( '\nFinished reading files...' )

    croppedGLIMap = gliMap[ :, topRow : bottomRow, leftCol : rightCol ]
    if saveImg == True: 
        if VERBOSE: print('\nSaving gli subregion map...')
        tiffy.imsave( outputDir + sampleName + '-gliMap.tif',croppedGLIMap.astype('uint16') )

    return croppedGLIMap

