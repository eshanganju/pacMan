"""Project module: abstracts 3D particles as 2D projections for CNN training

Authors:
	Eshan Ganju and Daniel Sinclair
	
"""

import numpy as np
import math
import tifffile as tf
from numba import jit

TESTING = False

def getTiffStacks2(labelledMap,saveData=True,fileName='',outputDir='',
					size=250,normalize=True):
	"""
	"""
	numParticles = int( labelledMap.max() )

	particleProjectionStack = np.zeros((1,size,size,3))

	for label in range( 1, numParticles + 1 ):
		print('Checking label:' + str(label) + '/' + str(numParticles))
		
		# extract volume around particle
		paddedNormParticle = cropAndPadParticle( labelledMap,label,size,
													saveData=True,
													fileName=fileName,
													outputDir=outputDir)

		if TESTING == True:
			saveFileName = outputDir + fileName + '-' + str(np.round(label)) +'-Image.tif'
			tf.imwrite( saveFileName,paddedNormParticle.astype('float') )
		
		# get projections
		particleProjectionStack = takeProjections( paddedNormParticle, 
																	size=size,
																	normalize=True )
		
		saveNameCurrent=outputDir + fileName + '-' + str(np.round(label)) +'-ImageStack.tif'
		
		tf.imwrite(saveNameCurrent,particleProjectionStack.astype('float'))



@jit(nopython=True)
def cropAndPadParticle(labelledMap,label,size, saveData=True,fileName='',outputDir=''):
	"""
	"""
	loc = np.where(labelledMap == label)
	dimZ = int(loc[0].max() - loc[0].min())
	dimY = int(loc[1].max() - loc[1].min())
	dimX = int(loc[2].max() - loc[2].min())
		
	padZ = int(np.round( ( size - dimZ ) / 2 ))
	padY = int(np.round( ( size - dimY ) / 2 ))
	padX = int(np.round( ( size - dimX ) / 2 ))
	

	croppedParticle = labelledMap[ loc[0].min() : loc[0].max()+1,
									loc[1].min() : loc[1].max()+1,
									loc[2].min() : loc[2].max()+1 ]

	paddedParticle = np.zeros( ( size,size,size ) )
	
	paddedParticle[ padZ : dimZ+padZ+1, padY : dimY+padY+1, padX : dimX+padX+1] = croppedParticle

	cleanedPaddedParticle = removeOtherLabels(paddedParticle, label)

	# convert label number to ones
	normalizedCleanedPaddedParticle = cleanedPaddedParticle//label

	return normalizedCleanedPaddedParticle

@jit(nopython=True)
def removeOtherLabels(paddedParticle, label):
	"""
	"""
	return np.where(np.logical_or(paddedParticle > label, paddedParticle < label), np.zeros_like(paddedParticle), paddedParticle)


@jit(nopython=True)
def takeProjections(paddedNormParticle,size,normalize):
	"""
	"""
	particleProjectionParticle = np.zeros((1,size,size,3))

	for direction in range(0,3):

		proj = np.sum( paddedNormParticle, axis=direction)
		
		if normalize == True:
			proj = proj/proj.max()
		
		particleProjectionParticle[:,:,:,direction] = proj


	return particleProjectionParticle
