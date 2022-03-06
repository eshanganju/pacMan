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

def extractSubregionFromTiffFile( fileLoc, Z, Y, X, lngt, calib, zReference = 'middle',xyReference='topLeft',invImg=False, saveImg=False,outputDir='', sampleName='' ):
	"""This function reads data from a single tiff stack
	The entire image is read and a cubic subregion specified by the user is returned

	Parameter
	----------
	fileLoc : str
		location of the image stack
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

	# zReference options
	if zReference == 'center':
		upperSlice = int(Z) + int( ( lngt / 2 ) / calib )
		lowerSlice = int(Z) - int( ( lngt / 2 ) / calib )

	elif zReference == 'low':
		upperSlice = int(Z) + int( lngt / calib ) 
		lowerSlice = int(Z) 

	elif zReference == 'high':
		upperSlice = int(Z)
		lowerSlice = int(Z) - int( lngt / calib )

	# xyRefence options
	if xyReference == 'topLeft':
		topRow = int(Y)
		bottomRow = int(Y) + int(lngt/calib)
		leftCol = int(X)
		rightCol = int(X) + int(lngt/calib)

	elif xyReference == 'center':
		topRow = int(Y) - int(lngt/calib)//2
		bottomRow = int(Y) + int(lngt/calib)//2
		leftCol = int(X) - int(lngt/calib)//2
		rightCol = int(X) + int(lngt/calib)//2

	gliMap = tiffy.imread( fileLoc )

	if invImg == True: gliMap = invertImage(gliMap)
	print( '\nFinished reading files...' )

	croppedGLIMap = gliMap[ lowerSlice : upperSlice + 1, topRow : bottomRow + 1, leftCol : rightCol + 1 ]

	if saveImg == True: 
		if VERBOSE: print('\nSaving gli subregion map...')
		tiffy.imsave( outputDir + sampleName + '-subMap.tif',croppedGLIMap.astype('uint16') )

	return croppedGLIMap

def _readTiffStack(fileLoc='',invImg=False,downSample=False,reduceBitDepth=False,initalBitDepth=0,finalBitDepth=0,saveImg=False,outputDir='',sampleName=''):
	"""
	"""
	image = tf.imread(fileLoc)

	# Down sample first
	if downSample == True: image = _downSampleData(initalImageData=image,saveImg=False)

	# reduce Bit Depth later
	if reduceBitDepth == True: image = _reduceBitDepthData(image,saveImg=False)

	# Save image - images will not be saved in the downSample and reduceBitDepth parts of the code
	if saveImg==True:
		print('... saving file.')

		if downSample == True:

			if reduceBitDepth == True:
				saveFileName = sampleName+'-tiffStack-BitDepthReduced-DownSampled.tif'
				tiffy.imwrite(saveFileName,image)

			else:
				saveFileName = sampleName+'-tiffStack-BitDepthReduced.tif'
				tiffy.imwrite(saveFileName,image)
		
		else:
			saveFileName = sampleName+'-tiffStack.tif'
			tiffy.imwrite(saveFileName,image)

		print('...file saved')

	else: print('... No file saved')

	return image



def _downSampleData(initalImageData,saveImg=False,outputDir='',sampleName=''):
	"""
	"""

def _reduceBitDepthData(intialImage,initalBitDepth=0,finalBitDepth=0,saveImg=False,outputDir='',sampleName=''):
	"""
	"""


def invertImage( gliMapToInvert ):
	"""Inverts the image  
	"""
	invertedGLIMap = np.flip(gliMapToInvert,axis=1)
	return invertedGLIMap
