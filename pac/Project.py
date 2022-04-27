"""Project module: abstracts 3D particles as 2D projections for CNN training

Authors:
	Daniel Sinclair
	Eshan Ganju
"""

import numpy as np
import math
import tifffile as tf


def getTiffstacks( labelledMap, saveData=False, fileName='',outputDir='', crop=True, size=250, normalize=True ):
	
	"""Creates 3 stacks of tiff images with each numbered slice containing a projection
	   of the corresponding labelled particle

	Parameters
	----------
	labelledMap : unsigned integer narray

	saveData : bool

	fileName : string

	outputDir : string

	crop : bool (recommended to reduce file size and improve efficiency)

	size : int

	normalize : bool

	Return:
	zProjections, yProjections, xProjections : unsigned float ndarrays
		Arrays with the x and y dimensions of labelledMap sliced along orthogonal axes
		Dimensions reduced to a square of size x size
		Each numbered slice corresponds to a labelled particle and contains its projection along Z, Y, or X axis
	"""
	
	numberOfParticles = 5#int(labelledMap.max())
	labelledShape = labelledMap.shape
	if crop == True:
		zProjections = np.zeros((numberOfParticles,size,size))
		yProjections = np.zeros((numberOfParticles,size,size))
		xProjections = np.zeros((numberOfParticles,size,size))
	else:
		zProjections = np.zeros((numberOfParticles,labelledShape[1],labelledShape[2]))
		yProjections = np.zeros((numberOfParticles,labelledShape[0],labelledShape[2]))
		xProjections = np.zeros((numberOfParticles,labelledShape[0],labelledShape[1]))

	for label in range(1,numberOfParticles+1):
		print('Computing projections for '+str(label)+'/'+str(numberOfParticles))
		zProjections[label-1,:,:]=getProjection(labelledMap, label, direction = 0, crop=crop, size=size, normalize=normalize)
		yProjections[label-1,:,:]=getProjection(labelledMap, label, direction = 1, crop=crop, size=size, normalize=normalize)
		xProjections[label-1,:,:]=getProjection(labelledMap, label, direction = 2, crop=crop, size=size, normalize=normalize)
		

	if saveData == True:
		print('Saving projections...')
		file = outputDir+fileName
		tf.imwrite(file+'_Zproject',zProjections)
		tf.imwrite(file+'_Yproject',yProjections)
		tf.imwrite(file+'_Xproject',xProjections)

	return zProjections, yProjections, xProjections

def getProjection( labelledMap, label, direction = 0, crop=True, size=250, normalize=True):

	"""Reduces the dimensions of a label in labelledMap from 3 to 2 using Numpy's reduce function

	Parameters
	----------
	labelledMap : unsigned integer array

	label : integer

	direction : integer (0-2)

	crop : bool

	size : integer

	normalize : bool

	Return
	projection : unsigned float ndarray
		2D array computed from a single label summed along the axis selected using direction

	"""
	if crop == True:
		particle = cropParticle(labelledMap, label, size)
	else:
		particle = np.zeros_like(labelledMap)
		particle[np.where(labelledMap ==label)] = 1

	
	projection = np.add.reduce(particle, axis=direction)
		
	if normalize == True:
		projection = projection/projection.max()

	return projection

def cropParticle(labelledMap, label, size):
	
	"""Isolates a label in a cube with dimensions equal to size.
	
	Parameters
	----------
	labelledMap : unsigned integer array

	label : integer

	size : integer

	Return
	paddedParticle : unsigned integer array
	Ndarray with the selected label centered and value = 1. Dimensions are size x size x size

	"""
	
	loc = np.where(labelledMap == label)
	dimZ = loc[0].max() - loc[0].min()
	dimY = loc[1].max() - loc[1].min()
	dimX = loc[2].max() - loc[2].min()
	padZ = round((size - dimZ)/2)
	padY = round((size - dimY)/2)
	padX = round((size - dimX)/2)
	

	croppedParticle = np.copy(labelledMap[ loc[0].min() : loc[0].max()+1,\
							  loc[1].min() : loc[1].max()+1,\
							  loc[2].min() : loc[2].max()+1 ])
	

	paddedParticle = np.zeros((size,size,size))
	paddedParticle[padZ : dimZ+padZ+1, padY : dimY+padY+1, padX : dimX+padX+1] = croppedParticle
	paddeParticle=paddedParticle/label
	return paddedParticle