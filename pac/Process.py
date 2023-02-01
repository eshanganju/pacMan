""" Process.py is a helper function that converts images to desired state
"""

import numpy as np
import tifffile as tf
from PIL import Image
from numba import jit
import matplotlib.pyplot as plt

TESTING = True
VERBOSE = True

def convertImageBitDepth(originalBitDepth,targetBitDepth,
							saveData=True,fileName='',outputDir=''):
	"""
	"""

def resizeData(originalSize,targetSize=None, targetRatio=1,
				saveData=True,fileName='',outputDir=''):
	"""Simple resize of image while retaining the aspect ratio
	"""

#@jit(nopython=True)
def alignGrainBoundaryMap(referenceGrainBoundaryMap, grainBoundaryMap, padding=100,
							saveData=True, fileName='', outputDir=''):
	"""Takes two grain boundary maps and aligns them

	Alignment is done by sliding the grain boundary map on top of the reference boundary map
	At different positions, the corespondence between the boundaries is checked

	The position with the maximum correspondence is taken as the aligned position
	"""
	# Convert images to binary with max value of 1 and min value of 0
	maxGLI_REF = referenceGrainBoundaryMap.max()
	maxGLI_GB = grainBoundaryMap.max()

	referenceMap = referenceGrainBoundaryMap//maxGLI_REF
	gbMap = grainBoundaryMap//maxGLI_GB

	# Apply padding to reference image
	referenceMapSize0 = referenceMap.shape[0]
	referenceMapSize1 = referenceMap.shape[1]

	paddedReferenceMapSize0 = referenceMap.shape[0] + 2*padding
	paddedReferenceMapSize1 = referenceMap.shape[1] + 2*padding

	paddedReferenceMap = np.zeros( ( paddedReferenceMapSize0, paddedReferenceMapSize1 ) )
	paddedReferenceMap[ padding:padding+referenceMapSize0, padding:padding+referenceMapSize1 ] = referenceMap 

	# GB Map size - for padding later
	gbMapSize0 = gbMap.shape[0]
	gbMapSize1 = gbMap.shape[1]

	# Set current max overlap value = 0
	currentMaxOverlapPixelCount = 0

	# Max overlap image and multipled maps
	currentMaxOverlapMap = np.zeros_like(paddedReferenceMap)
	currentMaxMultipliedMap = np.zeros_like(paddedReferenceMap)

	# Execute for-loop to check over lap
	for pad0 in range(0, padding+1):
		
		for pad1 in range(0, padding+1):

			if VERBOSE == True:
				print('pad0: ' + str(pad0) + ', pad1: ' + str(pad1))

			# Apply offset padding to GB map
			currentPaddedMap = np.zeros_like(paddedReferenceMap)
			currentPaddedMap[ pad0:pad0+gbMapSize0, pad1:pad1+gbMapSize1 ] = gbMap	
			
			# Multiply padded GB map and padded refernce map
			currentMulipliedMap = currentPaddedMap * paddedReferenceMap

			# Take sum of multipled image to get current overlap pixel count
			totalOverlappingPixels = currentMulipliedMap.sum()

			# If current overlap value > max overlap value, update Max overlap image
			if totalOverlappingPixels > currentMaxOverlapPixelCount:
				print('\tMax-overlap map updated')
				currentMaxOverlapMap = currentPaddedMap
				currentMaxMultipliedMap = currentMulipliedMap
				currentMaxOverlapPixelCount = totalOverlappingPixels

	# Save max overlap Image and reference image
	if saveData == True: 
		fileNameOverLap = outputDir + fileName + '-Overlap.tif'
		fileNameReference = outputDir + fileName + '-Reference.tif'
		fileNameGB = outputDir + fileName + '-GB.tif'

		tf.imwrite( fileNameReference,(paddedReferenceMap*maxGLI_REF).astype( 'uint8' ) )
		tf.imwrite( fileNameGB, (currentMaxOverlapMap*maxGLI_REF).astype( 'uint8' )  )
		tf.imwrite( fileNameOverLap,(currentMaxMultipliedMap*maxGLI_REF).astype( 'uint8' ) )

	return paddedReferenceMap, currentMaxOverlapMap, currentMaxMultipliedMap


