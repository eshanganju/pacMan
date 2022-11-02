"""Measure module: carries out measurements on the segmented data

Authors:
	Eshan Ganju
	Daniel Sinclair
"""

import numpy as np
import math
import statistics
import scipy
from numba import jit
from scipy.ndimage import binary_erosion as erode
import spam.label as slab
from uncertainties import unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *
from pac import Segment
from skimage.measure import marching_cubes, mesh_surface_area, regionprops
import tifffile as tf

VERBOSE = False				# Show all the text when running the functions
TESTING = True				# Keep if TESTING == True: Print('') while testing in function with "_" prefix 


def getPSDAll( labelledMap , calibrationFactor=1.0, getEqsp=True, getCaMax=True,
				getCaMed=True, getCaMin=True, getFeretMax=True, getFeretMin=True,):			
	"""This module returns the paricle size distribution for the 6 size
	parameters commonly used to quantify particle size.

	The reason this module exists is because different size parameters capture
	the size of the different sands more accurately. The parameter that best
	captures the size distribution obtained from sieve analysis should generally
	be used. The appropriate parameter depends on the shape of the particle.

	As particles crush, the parameter that best captures the size of the
	particles may change. Generally the "minimum feret diameter" works quite well
	in all cases.

	Parameters
	----------
	labelledMap : unsigned integer ndarray

	calibrationFactor : float

	getEqsp : bool

	getCaMax : bool

	getCaMed : bool

	getCaMin : bool

	getFeretMax : bool

	getFeretMin : bool

	saveData : bool

	sampleName : string

	outputDir : string

	Return
	------
	psdEqsp : floats n by 2 array

	psdCaMax : floats n by 2 array

	psdCaMed : floats n by 2 array

	psdCaMin : floats n by 2 array

	psdFeretMax : floats n by 2 array

	psdFeretMin : floats n by 2 array

	"""

	psArray = getParticleSizeArray( labelledMap,
									calibrationFactor=calibrationFactor,
									saveData=saveData,
									sampleName=sampleName,
									outputDir=outPutDir)

	psdEqsp = getParticleSizeDistribution( psArray, sizeParam='eqsp',
											sampleName=sampleName,
											saveData=saveData,
											outputDir=outputDir )

	psdCaMax = getParticleSizeDistribution( psArray, sizeParam='caMax',
											sampleName=sampleName,
											saveData=saveData,
											outputDir=outputDir )

	psdCaMed = getParticleSizeDistribution( psArray, sizeParam='caMed',
											sampleName=sampleName,
											saveData=saveData,
											outputDir=outputDir )

	psdCaMin = getParticleSizeDistribution( psArray, sizeParam='caMin',
											sampleName=sampleName,
											saveData=saveData,
											outputDir=outputDir )

	psdFeretMax = getParticleSizeDistribution( psArray, sizeParam='feretMax',
												sampleName=sampleName,
												saveData=saveData,
												outputDir=outputDir )

	psdFeretMin = getParticleSizeDistribution( psArray, sizeParam='feretMin',
												sampleName=sampleName,
												saveData=saveData,
												outputDir=outputDir )

	return psdEqsp, psdCaMax, psdCaMed, psdCaMin, psdFeretMax, psdFeretMin


def getParticleSizeArray( labelledMapForParticleSizeAnalysis, calibrationFactor=1,
							getCaDia=True,getFeretDia=True,
							saveData=False, sampleName='', outputDir=''):
	"""Computes particle size parameters for all the labels in the segmented data

	Parameters
	----------
	labelledMapForParticleSizeAnalysis : unsigned integer ndarray

	calibrationFactor : float

	saveData : bool

	sampleName : string

	outputDir : string

	Return:
	particleSizeDataSummary : unsigned float ndarray
		Particle size array containing columns
		[0] Label index, [1] Volume, [2] Eqsp, [3] Centroidal - max,
		[4] Centroidal - med, [5] Centroidal - min, [6] Feret-max X, [7] Feret -min X
	"""
	numberOfParticles = int( labelledMapForParticleSizeAnalysis.max() )

	particleSizeDataSummary = np.zeros( ( numberOfParticles + 1 , 8 ) )
	print( '\nStarting measurement of particles' )
	print( '------------------------------------*' )

	for particleNum in range( 1, numberOfParticles + 1 ):
		print( "Computing size of", particleNum, "/", numberOfParticles, "particle" )
		particleSizeDataSummary[particleNum, 0] = particleNum

		# Equivalent sphere diameters
		vol, eqspDia = getEqspDia( labelMap=labelledMapForParticleSizeAnalysis, label=int(particleNum) )
		particleSizeDataSummary[particleNum, 1] = vol * (calibrationFactor**3)
		particleSizeDataSummary[particleNum, 2] = eqspDia * calibrationFactor

		# Centroidal axes lengths
		if getCaDia == True:
			caMax, caMed, caMin = getPrincipalAxesLengths( labelMap=labelledMapForParticleSizeAnalysis,label=int(particleNum) )
			particleSizeDataSummary[particleNum, 3] = caMax * calibrationFactor
			particleSizeDataSummary[particleNum, 4] = caMed * calibrationFactor
			particleSizeDataSummary[particleNum, 5] = caMin * calibrationFactor

		# TODO: Fix the feret diameter issue
		"""
		The SPAM code runs and gives error in the apply transformation part of the code
		"""

		# Feret diameters
		# feretMax, feretMin = getMinMaxFeretDia( labelMap=labelled *-MapForParticleSizeAnalysis,
		#                                         label=int(particleNum), numOrts=100, numRots=10)

		if getFeretDia == True:
			feretDia = getFeretDiametersSPAM( lab=labelledMapForParticleSizeAnalysis, labelList=int(particleNum)  )
			feretMax = feretDia [0,0]
			feretMin = feretDia [0,1]
			particleSizeDataSummary[particleNum, 6] = feretMax * calibrationFactor
			particleSizeDataSummary[particleNum, 7] = feretMin * calibrationFactor

	# Removing the extra zero in the summary (first row)
	particleSizeDataSummary = np.delete( particleSizeDataSummary, 0, 0 )

	if saveData == True:
		if VERBOSE:  print('\nSaving particle size list...')
		np.savetxt( outputDir + sampleName + '-particleSizeList.csv',particleSizeDataSummary, delimiter=',')

	# [ Label, Volume(vx), Size0(px or mm), Size1(px or mm), Size2(px or mm), Size3(px or mm), Size4(px or mm), Size5(px or mm)]
	return particleSizeDataSummary


def getEqspDia( labelMap, label ):
	"""Calculates the equivalent sphere diameter of a particle in px units

	The equivalent sphere diameter is the diameter of the sphere with
	the same volume as the partcle

	Parameters
	----------
	labelMap : unsigned integer ndarray

	label : unsigned integer

	Return
	------
	volume : integer

	eqspLabelDia : float

	"""
	labOnlyMap = np.zeros_like( labelMap )
	labOnlyMap[np.where(labelMap == label)] = 1
	volume = labOnlyMap.sum()
	eqspLabelDia = (6*volume/(math.pi))**(1/3)

	return volume, eqspLabelDia


def getPrincipalAxesOrtTable( labelMapForParticleOrientation, saveData=False, sampleName='', outputDir=''):
	"""Computes a table of the major axes of the particles - orientation of the particles long axis
	"""
	numberOfParticles = int( labelMapForParticleOrientation.max() )
	particleOrientationDataSummary = np.zeros( ( numberOfParticles + 1 , 4 ) )  # This is needed because of ZYX + the particle number
	
	for particleNum in range( 1, numberOfParticles + 1 ):
		print( "Computing orientation of", particleNum, "/", numberOfParticles, "particle" )
		particleOrientationDataSummary[particleNum, 0] = particleNum
		ortZZ, ortYY, ortXX = getMajorPrincipalAxesOrt( labelMap=labelMapForParticleOrientation, label=int(particleNum) )
		particleOrientationDataSummary[particleNum, 1] = ortZZ
		particleOrientationDataSummary[particleNum, 2] = ortYY
		particleOrientationDataSummary[particleNum, 3] = ortXX

	# Removing the extra zero in the summary (first row)
	particleOrientationDataSummary = np.delete( particleOrientationDataSummary, 0, 0 )

	# Saving files
	if saveData == True:
		if VERBOSE:  print('\nSaving particle ort list...')
		np.savetxt( outputDir + sampleName + '-particleOrtList.csv',particleOrientationDataSummary, delimiter=',')

	# [ Label, Volume(vx), Size0(px or mm), Size1(px or mm), Size2(px or mm), Size3(px or mm), Size4(px or mm), Size5(px or mm)]
	return particleOrientationDataSummary


def getMajorPrincipalAxesOrt( labelMap, label ):
	"""Computes the orientation of the principal axes along which the particle has the "largest length"
	The eigen vector corresponding to the largest eigenvalue is the vector along which the point cloud of the particle has the largest length.
	"""
	zyxofLabel = getZYXLocationOfLabel(labelMap,label)

	covarianceMatrix = np.cov( zyxofLabel.T )

	eigval, eigvec = np.linalg.eig( covarianceMatrix )

	# Eigenvector corresponding to max eigenvalue
	maxEigvalEigvecIndex = int(np.argmax(eigval))

	ortZZ = eigvec[0,maxEigvalEigvecIndex]
	ortYY = eigvec[1,maxEigvalEigvecIndex]
	ortXX = eigvec[2,maxEigvalEigvecIndex]

	return ortZZ, ortYY, ortXX


def getPrincipalAxesLengths( labelMap, label ):
	"""Computes the principal axes lengths of the particle in px units
	
	Parameters
	----------
	labelMap : unsigned integer ndarray
	
	label : unsigned integer
	
	Return
	------
	centroidalAxesLengthMax : unsigned float
	
	centroidalAxesLengthMed : unsigned float
	
	centroidalAxesLengthMin : unsigned float
	"""
	zyxofLabel = getZYXLocationOfLabel(labelMap,label)
	
	covarianceMatrix = np.cov( zyxofLabel.T )
	#print('\tobtained COV matrix')
	eigval, eigvec = np.linalg.eig( covarianceMatrix )
	#print('\tobtained eigvectors matrix')
	meanZ, meanY, meanX = getCenterOfGravityFromZYXLocations( zyxLocationData=zyxofLabel )
	
	meanMatrix = np.zeros_like( zyxofLabel )
	
	meanMatrix[ :, 0 ] = meanZ
	meanMatrix[ :, 1 ] = meanY
	meanMatrix[ :, 2 ] = meanX
	
	centeredLocationData = zyxofLabel - meanMatrix
	
	rotationMatrix = np.zeros( ( 3, 3 ) )
	
	rotationMatrix[ :, 0 ] = eigvec[ 0 ]
	rotationMatrix[ :, 1 ] = eigvec[ 1 ]
	rotationMatrix[ :, 2 ] = eigvec[ 2 ]
	
	rotCentPointCloud = ( np.matmul( rotationMatrix, centeredLocationData.T ) ).T
	
	caDims = np.zeros( ( 3, 1 ) )
	caDims[ 0 ] = rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 0 ].min()
	caDims[ 1 ] = rotCentPointCloud[ :, 1 ].max() - rotCentPointCloud[ :, 1 ].min()
	caDims[ 2 ] = rotCentPointCloud[ :, 2 ].max() - rotCentPointCloud[ :, 2 ].min()
	
	centroidalAxesLengthMax = max( caDims )[ 0 ]
	centroidalAxesLengthMin = min( caDims )[ 0 ]
	centroidalAxesLengthMed = statistics.median( caDims )[ 0 ]
	return centroidalAxesLengthMax, centroidalAxesLengthMed, centroidalAxesLengthMin


# TODO: This is messy - needs to be cleaned
def getFeretDiametersSPAM(lab, labelList=None, boundingBoxes=None, centresOfMass=None, numberOfOrientations=100,
							margin=0, interpolationOrder=0, returnOrts = False):
	"""
	Calculates (binary) feret diameters (caliper lengths) over a number of equally-spaced orientations
	and returns the maximum and minimum values, as well as the orientation they were found in.

	Parameters
	----------
		lab : 3D array of integers
			Labelled volume, with lab.max() labels

		labelList: list of ints, optional
			List of labels for which to calculate feret diameters and orientations. Labels not in lab are ignored. Outputs are given in order of labelList.
			If not defined (Default = None), a list is created from label 0 to lab.max()

		boundingBoxes : lab.max()x6 array of ints, optional
			Bounding boxes in format returned by ``boundingBoxes``.
			If not defined (Default = None), it is recomputed by running ``boundingBoxes``

		centresOfMass : lab.max()x3 array of floats, optional
			Centres of mass in format returned by ``centresOfMass``.
			If not defined (Default = None), it is recomputed by running ``centresOfMass``

		numberOfOrientations : int, optional
			Number of trial orientations in 3D to measure the caliper lengths in.
			These are defined with a Saff and Kuijlaars Spiral.
			Default = 100

		margin : int, optional
			Number of pixels by which to pad the bounding box length to apply as the margin in spam.label.getLabel().
			Default = 0

		interpolationOrder = int, optional
			Interpolation order for rotating the object.
			Default = 0

	Returns
	-------
		feretDiameters : lab.max()x2 (or len(labelList)x2 if labelList is not None) array of integers
			The max and min values of the caliper lengths of each labelled shape.
			Expected accuracy is +- 1 pixel

		feretOrientations : lab.max()x6 (or len(labelList)x6 if labelList is not None) array of floats
			2 x Z,Y,X components of orientations of the max and min caliper lengths

	Notes
	-----
		Function contributed by Estefan Garcia (Caltech, previously at Berkeley)
	"""

	#Notes
	#-------
		#Must import spam.DIC to use this function because it utilizes the computePhi and applyPhi functions.
		#This function currently runs in serial but can be improved to run in parallel.

	import spam.DIC
	import spam
	import spam.label as slab
	import spam.plotting as splt
	import spam.deformation as sdef

	#print(1)

	labelType = '<u4'
	lab = lab.astype(labelType)

	if labelList is None:
		labelList = list(range(0,lab.max()+1))
		feretDiameters = np.zeros((lab.max() + 1, 2))
		feretOrientations = np.zeros((lab.max()+1,6))

	elif type(labelList) is not list and type(labelList) is not np.ndarray:
		# Allow inputs to be ints or of type np.ndarray
		labelList = [labelList]
		feretDiameters = np.zeros((len(labelList),2))
		feretOrientations = np.zeros((len(labelList),6))

	else:
		feretDiameters = np.zeros((len(labelList),2))
		feretOrientations = np.zeros((len(labelList),6))

	#print('Calculating Feret diameters for '+str(len(labelList))+' label(s).')

	if boundingBoxes is None:
		boundingBoxes = spam.label.boundingBoxes(lab)

	if centresOfMass is None:
		centresOfMass = spam.label.centresOfMass(lab, boundingBoxes=boundingBoxes)

	# Define test orientations
	testOrientations = splt.orientationPlotter.SaffAndKuijlaarsSpiral(4*numberOfOrientations)

	i=0
	while i < len(testOrientations):
		if (testOrientations[i] < 0).any():
			testOrientations = np.delete(testOrientations,i,axis=0)
		else:
			i+=1

	#print(2)

	# Compute rotation of trial orientations onto z-axis
	rot_axes = np.cross(testOrientations,[1.,0.,0.])
	rot_axes/=np.linalg.norm(rot_axes,axis=1,keepdims=True)
	theta=np.reshape(np.rad2deg(np.arccos(np.dot(testOrientations,[1.,0.,0.]))),[len(testOrientations),1])

	# Compute Phi and its inverse for all trial orientations
	Phi = np.zeros((len(testOrientations),4,4))
	transf_R = rot_axes*theta

	for r in range(0,len(transf_R)):
		transformation = {'r': transf_R[r]}
		Phi[r] = sdef.deformationFunction.computePhi(transformation)

	Phi_inv = np.linalg.inv(Phi)

	#print(3)

	# Loop through all labels provided in labelList. Note that labels might not be in order.
	for labelIndex in range(0,len(labelList)):
		#print(4)

		label = labelList[labelIndex]

		if label in lab and label > 0: #skip if label does not exist or if zero

			particle = slab.label.getLabel(lab,
											label,
											boundingBoxes= boundingBoxes,
											centresOfMass= centresOfMass,
											extractCube= True,
											margin= margin,
											maskOtherLabels= True)

			subvol = particle['subvol']

			# Initialize DMin and DMax using the untransformed orientation
			subvol_transformed_BB = spam.label.boundingBoxes(subvol > 0.5)
			zWidth = subvol_transformed_BB[1,1] - subvol_transformed_BB[1,0] + 1
			yWidth = subvol_transformed_BB[1,3] - subvol_transformed_BB[1,2] + 1
			xWidth = subvol_transformed_BB[1,5] - subvol_transformed_BB[1,4] + 1

			index_max = np.argmax([zWidth,yWidth,xWidth])
			index_min = np.argmin([zWidth,yWidth,xWidth])

			DMax = max([zWidth, yWidth, xWidth])
			DMin = min([zWidth, yWidth, xWidth])

			maxOrientation = [np.array([1.,0.,0.]),
								np.array([0.,1.,0.]),
								np.array([0.,0.,1.])][index_max]

			minOrientation = [np.array([1.,0.,0.]),
								np.array([0.,1.,0.]),
								np.array([0.,0.,1.])][index_min]

			#print(5)

			for orientationIndex in range(0,len(testOrientations)):
				# Apply rotation matrix about centre of mass of particle
				subvol_centreOfMass = spam.label.centresOfMass(subvol)

				#print('\t Ort ' + str(orientationIndex) + ' 5-1')
				subvol_transformed = spam.DIC.applyPhi(subvol,
														Phi = Phi[orientationIndex],
														PhiCentre = subvol_centreOfMass[1],
														interpolationOrder=interpolationOrder)

				#print('\t Ort ' + str(orientationIndex) + ' 5-2')

				# Use bounding box of transformed subvolume to calculate particle widths in 3 directions
				subvol_transformed_BB = spam.label.boundingBoxes(subvol_transformed > 0.5)
				zWidth = subvol_transformed_BB[1,1] - subvol_transformed_BB[1,0] + 1
				yWidth = subvol_transformed_BB[1,3] - subvol_transformed_BB[1,2] + 1
				xWidth = subvol_transformed_BB[1,5] - subvol_transformed_BB[1,4] + 1

				#print('\t Ort ' + str(orientationIndex) + ' 5-3')

				# Check if higher than previous DMax or lower than previous DMin
				index_max = np.argmax([DMax,zWidth,yWidth,xWidth])
				index_min = np.argmin([DMin,zWidth,yWidth,xWidth])
				DMax = max([DMax,zWidth,yWidth,xWidth])
				DMin = min([DMin,zWidth,yWidth,xWidth])

				# Update orientations for DMax and DMin
				maxOrientation = [maxOrientation,
									testOrientations[orientationIndex],
									np.matmul(Phi_inv[orientationIndex,:3,:3],np.array([0,1,0])),
									np.matmul(Phi_inv[orientationIndex,:3,:3],np.array([0,0,1]))][index_max]

				minOrientation = [minOrientation,
									testOrientations[orientationIndex],
									np.matmul(Phi_inv[orientationIndex,:3,:3],np.array([0,1,0])),
									np.matmul(Phi_inv[orientationIndex,:3,:3],np.array([0,0,1]))][index_min]

				#print(6)

			feretDiameters[labelIndex,:] = [DMax,DMin]
			feretOrientations[labelIndex,:] = np.concatenate([maxOrientation,minOrientation])

	#print(7)

	if returnOrts == True: return feretDiameters,feretOrientations
	if returnOrts == False: return feretDiameters


def getMinMaxFeretDia( labelMap, label, numOrts=100, numRots=10 ):
	"""Computes the minimum and the maximum feret diameters of a particle
	
	The voxels of the particle are rotated along different unit vectors,
	which are obtained from the method proposed by Saff and Kuijlaar.
	After rotation, the lengths of the particle along the 3 axes are
	measured and the minimum and the maximum lengths measured are returned
	
	Parameters
	----------
	labelMap : unsigned integer ndarray
	
	label : unsigned integer
	
	numOrts: unsigned integer
		number of orientations about which particle will be rotated
	
	numRots: unsigned integer
		number of rotation angles about each orientation
	
	Return
	------
	minFeretDia : float
	
	maxFeretDia : float
	
	NOTES
	-----
	
	This code giced smaller particle sizes than the SPAM code. All analysis is reverted back to SPAM code till this issue is fixed. 
	"""
	if TESTING: print('\nChecking feret diameters of ' + str(np.round(label)) + '/' + str( np.round( labelMap.max() ) ) )
	
	zyxLocations = getZYXLocationOfLabel(labelMap, label).astype('float')
	
	cogZ, cogY, cogX = getCenterOfGravityFromZYXLocations ( zyxLocationData=zyxLocations )
	
	surfaceOnlyLabMap = getSurfaceVoxelsofParticle(labelMap, label, returnWithLabel = True)
	
	# zyxSurfaceLocations = getZYXLocationOfLabel(surfaceOnlyLabMap, label)
	zyxSurfaceLocations = zyxLocations
	print(zyxSurfaceLocations)
	
	meanMatrix = np.zeros_like(zyxSurfaceLocations)
	meanMatrix[:,0] = cogZ
	meanMatrix[:,1] = cogY
	meanMatrix[:,2] = cogX
	centeredZYXSurfaceLocations = zyxSurfaceLocations - meanMatrix
	
	unitVectors = getOrientationUsingSandK( numberOfOrientations=numOrts)
	anglesInRad = np.arange( 0, math.pi+0.00001, math.pi/numRots )
	
	if TESTING:
		print('Orientaions in radians:')
		print(anglesInRad)
	
	minFeretDiameter = 0
	maxFeretDiameter = 0
	
	for i in range(0,unitVectors.shape[0]):
		if TESTING: print('\tChecking along orientation: ' + str(unitVectors[i]))
		a = unitVectors[ i, 0 ]
		b = unitVectors[ i, 1 ]
		c = unitVectors[ i, 2 ]
		d = ( b**2 + c**2 ) ** 0.5

		Rz = np.identity(3)
		Rz[ 1, 1 ] = c/d
		Rz[ 1, 2 ] = -b/d
		Rz[ 2, 1 ] = b/d
		Rz[ 2, 2 ] = c/d

		RzInv = np.identity(3)
		RzInv[ 1, 1 ] = c/d
		RzInv[ 1, 2 ] = b/d
		RzInv[ 2, 1 ] = -b/d
		RzInv[ 2, 2 ] = c/d

		Ry = np.identity(3)
		Ry[ 0, 0 ] = d
		Ry[ 0, 2 ] = -a
		Ry[ 2, 0 ] = a
		Ry[ 2, 2 ] = d

		RyInv = np.identity(3)
		RyInv[ 0, 0 ] = d
		RyInv[ 0, 2 ] = a
		RyInv[ 2, 0 ] = d
		RyInv[ 2, 2 ] = -a

		preRotatedZYXSurfaceLocations = np.matmul( Ry, np.matmul( Rz, centeredZYXSurfaceLocations.T ) )

		for angle in anglesInRad:
			if TESTING: print('\t\tChecking at rotation of ' + str(angle) + ' radians')
			Rx = np.identity( 3 )
			Rx[ 0, 0 ] = math.cos( angle )
			Rx[ 0, 1 ] = -math.sin( angle )
			Rx[ 1, 0 ] = math.sin( angle )
			Rx[ 1, 1 ] = math.cos( angle )

			rotatedZYXSurfaceLocations = np.matmul( Rx, preRotatedZYXSurfaceLocations )
			correctedRotatedZYXSurfaceLocations = np.matmul( RzInv,np.matmul( RyInv, rotatedZYXSurfaceLocations ) )

			sizeRange = np.zeros( ( 3, 1 ) )
			# The one pixel value is added to account for the length of the pixel not considered
			sizeRange[ 0, 0 ] = correctedRotatedZYXSurfaceLocations[ 0,: ].max() - correctedRotatedZYXSurfaceLocations[ 0, : ].min() + 1
			sizeRange[ 1, 0 ] = correctedRotatedZYXSurfaceLocations[ 1,: ].max() - correctedRotatedZYXSurfaceLocations[ 1, : ].min() + 1
			sizeRange[ 2, 0 ] = correctedRotatedZYXSurfaceLocations[ 2,: ].max() - correctedRotatedZYXSurfaceLocations[ 2, : ].min() + 1

			if minFeretDiameter != 0:
				if minFeretDiameter > sizeRange.min():
					minFeretDiameter = sizeRange.min()
					if TESTING: print('\t\t\tUpdated minimum feret diameter to: ' + str(np.round(minFeretDiameter)))

			if maxFeretDiameter != 0:
				if maxFeretDiameter < sizeRange.max():
					maxFeretDiameter = sizeRange.max()
					if TESTING: print('\t\t\tUpdated max feret diameter to: ' + str(np.round(maxFeretDiameter)))

			if minFeretDiameter == 0:
				minFeretDiameter = min( sizeRange )
				if TESTING: print('\t\t\tInitial minimum feret diameter: ' + str(np.round(minFeretDiameter)))

			if maxFeretDiameter == 0:
				maxFeretDiameter = max( sizeRange )
				if TESTING: print('\t\t\tInitial max feret diameter: ' + str(np.round(maxFeretDiameter)))

	return maxFeretDiameter,  minFeretDiameter


def getCenterOfGravityFromZYXLocations( zyxLocationData ):
	"""Calculates the center of gravity of the particle from the
	zyx locations of the particles voxels

	Since the particle is made up of voxels, each of which have
	the same volume, and assuming the same density, the same mass,
	the center of gravity along each axis is the average location
	of the voxels.

	Parameters
	----------
	zyxLocationData : unsigned integer ndarray

	Return
	------
	centerOfGravityZ : float
	centerOfGravityY : float
	centerOfGravityX : float
	"""
	# centering
	centerOfGravityZ = np.average( zyxLocationData[ :, 0 ] )
	centerOfGravityY = np.average( zyxLocationData[ :, 1 ] )
	centerOfGravityX = np.average( zyxLocationData[ :, 2 ] )

	return centerOfGravityZ, centerOfGravityY, centerOfGravityX


def getSurfaceVoxelsofParticle( labelMap, label, returnWithLabel = True):
	"""Gets the surface voxels of the particle

	Rotating all the voxels of the particles is computationally expensive
	To measure size of the particle, only the surface voxels need to be rotated

	This function extracts the surface voxels of the particles by eroding the
	particle and then subtracting the eroded particle from the original
	paritcle.

	Parameters
	----------
	labelMap : unsigned integer ndarray

	label : unsigned integer

	returnWithLabel : bool
		Set this to true if you want to return the particle with each voxel
		assigned its label. If False, it returns with each voxel assigned 1

	Returns
	-------
	surfaceMap : unsigned integer ndarray
		Only the surface voxels of the particle
	"""
	originalParticleMap = np.zeros_like( labelMap )
	originalParticleMap[np.where(labelMap == label)] = 1
	erodedParticleMap = erode(originalParticleMap)
	surfaceMap = originalParticleMap - erodedParticleMap

	if returnWithLabel == False: return surfaceMap
	if returnWithLabel == True: return surfaceMap*int(label)


def getOrientationUsingSandK(numberOfOrientations=100):
	"""Computes the unit vectors distributed on the surface
	of the a unit sphere

	This is computed using the approach proposed by Saff and Kuijlaars.
	This code is similar to the code in SPAM.

	Parameters
	----------
	numberOfOrientations : integer

	Return
	------
	points : float ndarray
		Contains in each row the unit vectors distributed around a unit sphere
	"""
	M = int(numberOfOrientations)*2
	s = 3.6 / math.sqrt(M)
	delta_z = 2 / float(M)
	z = 1 - delta_z/2

	longitude = 0

	points = np.zeros( (numberOfOrientations,3) )

	for k in range( numberOfOrientations ):
		r = math.sqrt( 1 - z*z )
		points[k,2] = math.cos( longitude ) * r     #X
		points[k,1] = math.sin( longitude ) * r     #Y
		points[k,0] = z                             #Z
		z = z - delta_z
		longitude   = longitude + s/r

	return points


def getParticleSizeDistribution( psSummary, sizeParam='feretMin', sampleName='', saveData=True, outputDir='' ):
	"""Generates the particle size distribution from list of labels and sizes

	Different size parametres can be used to generate the grain size distribution.
	The size that most accurately matches the size distribution from the sieve
	analysis should ideally be used for the assessment of particle size distribution

	Parameters
	----------
	psSummary : n by 8 np array
		This contains the results of the getParticleSize function. The array
		should have rows equal to the number of particles in the samples,
		and one column each for (0) Label, (1) Volume, (2) equivalent spere
		diameter, (3) maximum centroidal axes length, (4) intermediate
		centroidal axes length, (5) minimum centroidal axes length, (6)
		maximum feret diameter, and (7) minimum feret diameter.
	sizeParam : string
		'eqsp' - for equivalent sphere diameter
		'caMax' - for max centroidal axes length
		'caMed' - for intermediate centroidal axes length
		'caMin' - for minimum centroidal axes length
		'feretMax' - for max feret diameter
		'feretMin' - for min feret diameter
	sampleName : string
	saveData : bool
	outputDir : string

	Return
	------
	gsdPP : n by 2 np array
		x column is the particle size and y column is the percentage passing
	"""
	if VERBOSE: print( '\nGetting GSD for assuming ' + str( sizeParam ) )

	if sizeParam == 'eqsp' : sizeCol = int(2)
	elif sizeParam == 'caMax' : sizeCol = int(3)
	elif sizeParam == 'caMed' : sizeCol = int(4)
	elif sizeParam == 'caMin' : sizeCol = int(5)
	elif sizeParam == 'feretMax' : sizeCol = int(6)
	elif sizeParam == 'feretMin' : sizeCol = int(7)

	label = psSummary[ : , 0 ].reshape( psSummary.shape[ 0 ] , 1 )
	vol = psSummary[ : , 1 ].reshape( psSummary.shape[ 0 ] , 1 )
	size = psSummary[: , sizeCol ].reshape( psSummary.shape[ 0 ] , 1 )
	gss = np.append( label, vol , 1 ).reshape( psSummary.shape[ 0 ] , 2 )
	gss = np.append( gss , size , 1 ).reshape( psSummary.shape[ 0 ] , 3 )

	gss = gss[ np.argsort( gss[ : , 2 ] ) ]

	totalVol = np.sum( gss[ : , 1 ] )
	pp = ( np.cumsum( gss[ : , 1 ] ) / totalVol * 100 ).reshape( gss.shape[ 0 ] , 1 )

	gsdPP = np.append( gss , pp, 1 )

	if saveData == True:
		if VERBOSE: print('\nSaving particle size distribution...')
		np.savetxt( outputDir + sampleName + '-'+ sizeParam +'-particleSizeDist.csv', gsdPP, delimiter=',')

	return gsdPP


def computeVolumeOfLabel( labelMap, label ):
	"""Calculated the volume of a label by counting number of voxels in label

	Parameters
	----------
	labelMap : unsigned integer ndarray
	label : unsigned integer

	Return
	------
	volumeOfLabel : unsigned integer
		Total volume of the label in px units
	"""
	labelOnlyMap = np.zeros_like(labelMap)
	labelOnlyMap[np.where(labelMap == label)] = 1
	volumeOfLabel = labelOnlyMap.sum()
	return volumeOfLabel


def getZYXLocationOfLabel( labelMap, label ):
	"""Obtains the zyx locations of the voxels of a particle

	Parameters
	----------
	labelMap : unsigned integer ndarray

	label : unsigned integer

	Return
	------
	zyxLocationData : unsigned integer ndarray
	"""
	particleLocationArray = np.where(labelMap == label)
	zyxLocationData = np.zeros( ( particleLocationArray[ 0 ].shape[ 0 ], 3 ) )
	zyxLocationData[:,0] = particleLocationArray[0]
	zyxLocationData[:,1] = particleLocationArray[1]
	zyxLocationData[:,2] = particleLocationArray[2]
	return zyxLocationData


def getRelativeBreakageHardin( psdOriginal, psdCurrent, smallSizeLimit=0.075, 
								saveData=True, sampleName='', outputDir='' ):
	"""Computes the relative breakage parameter according to the defintion by Hardin(1985)

	Hardin proposes that after a particle reaches a certain threshold size
	(sand-silt boundary), it will not break any further

	This function computes a relative breakage parameter that follows
	Hardin's proposal. Relative breakage parameters assume that the sand has an
	inital and an ultimate particle size distribution. In its uncrushed state the sand
	is in its intial particle size distribution and after undergoing the maximum
	curshing possible, it reaches its ultimate particle size distribution. When the
	sand is at its inital PSD, the relative breakage parameter is 0 and when
	it is at its ultimate PSD, the relative breakage parameter is 1 (or 100%).

	The relative breakage parameters is computed as a ratio of the "Current
	Breakage" of the sand to the "Potential Breakage" of the sand. The current
	breakage is computed as the area between the current PSD and the original
	PSD. The potential breakage is computed as the area between the ultimate PSD
	and the original PSD.

	According to Hardin, the ultimate PSD is a vertical line at the sand-silt boundary.

	Parameters
	----------
	psdOriginal : n by 2 numpy array
		Original particle size distribution with particle size in col 0 and percentage
		passing in col 1. the largest particle size is controlled by this gradation

	psdCurrent : n by 2 numpy array
		Current particle size distribution with particle size in col 0 and percentage
		passing in col 1. the largest particle size is controlled by this gradation

	smallSizeLimit : unsigned float
		Lower limit of integration (mm)

	Return
	------
	potentialBreakage : unsigned float
		Area in mm2 between the ultimate and original gradation

	currentBreakage : unsigned float
		Area in mm2 between the current and original gradation

	relativeBreakage : unsigned float
		The relative breakage parameter according to Hardin (1985)
	"""
	largeSizeLimit = psdOriginal[ :,0 ].max()

	areaUnderOriginalPSD = getAreaUnderPSDCurve( psdOriginal,
													maxSize=largeSizeLimit )

	areaUnderCurrentPSD = getAreaUnderPSDCurve( psdCurrent,
													maxSize=largeSizeLimit )

	areaUnderUltimatePSD = ( math.log10( largeSizeLimit ) -math.log10( smallSizeLimit ) ) * psdOriginal[ :,1 ].max()

	potentialBreakage = areaUnderUltimatePSD - areaUnderOriginalPSD
	currentBreakage = areaUnderCurrentPSD - areaUnderOriginalPSD
	relativeBreakage = currentBreakage/potentialBreakage*100

	if saveData==True:
		fileNameToSave = outputDir + sampleName + '-hrdRelBreakParams.csv'
		dataToSave = np.array( [ potentialBreakage,
									currentBreakage,
									relativeBreakage ] ).reshape(3,1)

		np.savetxt( fileNameToSave, dataToSave, delimiter=',')

	return potentialBreakage, currentBreakage, relativeBreakage


def getRelativeBreakageEinav( psdOriginal, psdCurrent , fracDim=2.6, smallSizeLimit=0.001,
								saveData=True, sampleName='',outputDir='',returnUPSD=True):
	"""Computes the relative breakage parameter according to the defintion of Einav (2007)

	Einav proposes that the largest particle size does not change with crushing and
	that ultimately, the crushed sand follows a fractal gradation curve

	This function  computes a relative breakage parameter that follows
	Einav's proposal. Relative breakage parameters assume that the sand has an
	inital and an ultimate particle size distribution. In its uncrushed state the sand
	is in its intial particle size distribution and after undergoing the "maximum
	curshing possible," it reaches its ultimate particle size distribution. When the
	sand is at its inital PSD, the relative breakage parameter is 0 and when
	it is at its ultimate PSD, the relative breakage parameter is 1 (or 100%).

	The relative breakage parameters is computed as a ratio of the "Current
	Breakage" of the sand to the "Potential Breakage" of the sand. The current
	breakage is computed as the area between the current PSD and the original
	PSD. The potential breakage is computed as the area between the ultimate PSD
	and the original PSD.

	According to Einav, the ultimate PSD is a fractal gradation that has a fractal
	dimension of 2.5-2.6

	Parameters
	----------
	psdOriginal : n by 2 numpy array
		Original particle size distribution with particle size in col 0 and percentage
		passing in col 1. the largest particle size is controlled by this gradation

	psdCurrent : n by 2 numpy array
		Current particle size distribution with particle size in col 0 and percentage
		passing in col 1.

	fracDim : unsigned float
		This is the fractal dimension used to compute the ultimate gradation

	smallSizeLimit : unsigned float
		Lower limit of integration (mm)

	Return
	------
	potentialBreakage : unsigned float
		Area in mm2 between the ultimate and original gradation

	currentBreakage : unsigned float
		Area in mm2 between the current and original gradation

	relativeBreakage : unsigned float
		The relative breakage parameter according to Einav (2007)
	"""
	largeSizeLimitOriginal = psdOriginal[ :,0 ].max()
	largeSizeLimitCurrent = psdCurrent[:,0].max()

	dmax = min(largeSizeLimitOriginal,largeSizeLimitCurrent)

	largeSizeLimit = psdOriginal[ :,0 ].max()

	areaUnderOriginalPSD = getAreaUnderPSDCurve( psdOriginal,
													maxSize=largeSizeLimit,
													minSize=smallSizeLimit,
													minOrZero='zero')

	areaUnderCurrentPSD = getAreaUnderPSDCurve( psdCurrent,
												maxSize=largeSizeLimit,
												minSize=smallSizeLimit,
												minOrZero='zero')

	psdUltimate = getUltimateFractalParticleSizeDistribution( minSize=smallSizeLimit,
																maxSize=dmax,
																fractalDimension=fracDim )

	areaUnderUltimatePSD = getAreaUnderPSDCurve( psdUltimate,
													maxSize=largeSizeLimit,
													minSize=smallSizeLimit,
													minOrZero='zero')

	potentialBreakage = areaUnderUltimatePSD - areaUnderOriginalPSD
	currentBreakage = areaUnderCurrentPSD - areaUnderOriginalPSD
	relativeBreakage = currentBreakage/potentialBreakage*100

	if saveData==True:
		fileNameToSave = outputDir + sampleName + '-envRelBreakParams.csv'
		dataToSave = np.array( [ potentialBreakage,
									currentBreakage,
									relativeBreakage ] ).reshape(3,1)

		np.savetxt( fileNameToSave, dataToSave, delimiter=',')

	return psdUltimate, potentialBreakage, currentBreakage, relativeBreakage


def getUltimateFractalParticleSizeDistribution(minSize=0.0,maxSize=0.0,fractalDimension=2.6,num=100):
	"""Gets the ultimate particle size distribution following the fractal
	distribution

	The equation is:
		P = (d/dMax)^(3-Xi)
			where p is the percentage by mass smaller than particle size d
			dMax is the largest particle size in the original sample
			Xi is the fractal dimension of the PSD in the ultimate state of crushing

	Parameters
	----------
	minSize : unsigned float
		smallest particle size in the PSD

	maxSize : unsigned float
		largest particle size in the PSD

	fractalDimension : unsigned float
		the fractal dimension of the ultimate gradation

	num : integer
		Number of points in the PSD

	Return
	------
	ultimatePSDFractal : unsigned float array
		Ultimate particle size distribution with particle size in col 0 and percentage
		passing in col 1. Col 0 is ascending order.
	"""
	if minSize == 0.0:
		minSize = float( input( 'Enter the min particle size (mm):' ) )
	if maxSize == 0.0:
		maxSize = float( input( 'Enter the max particle size (mm):' ) )

	ultimatePSDFractal = np.zeros((num,2))
	ultimatePSDFractal[ :, 0 ] = np.linspace( np.log10( minSize), np.log10( maxSize ), num )
	ultimatePSDFractal[ :, 0 ] = 10**ultimatePSDFractal[ :, 0 ]
	ultimatePSDFractal[ :, 1 ] = (( ultimatePSDFractal[ :, 0 ] / maxSize ) **( 3-fractalDimension ))*100

	return ultimatePSDFractal


def getAreaUnderPSDCurve( psd, maxSize=0.0, minSize=0.075, minOrZero='min'):
	"""Computes the area under the particle size distribution curve passed to it

	Since the GSD curves are plotted on a semi-log graph, the log10 of the x
	axis (particle size) is taken before computing the area.

	Parameters
	----------
	psd : array of floats
		Particle size distribution containing particle size (mm) in col 0 and
		percentage passing (%) in col 1. Col 0 should be in ascending order

	maxSize : unsigned float
		upper limit of integration in mm

	Return
	------
	areaUnderPSDCurve : unsigned float
		area under the psd
	"""
	if maxSize == 0.0: maxSize = float(input('Enter the max particle size (mm): '))
	psd = np.append( psd, np.array( [maxSize , 100.0] ).reshape( 1, 2 ), 0 )

	psdClipped = psd[np.where(psd[:,0] > minSize)]
	if minOrZero == 'min': minpp = min(psdClipped[:,1])
	if minOrZero == 'zero': minpp = 0
	val = np.array([minSize,minpp]).reshape(1,2)

	psd = np.concatenate((val,psdClipped), axis=0)

	numberOfPoints = psd.shape[0]
	areaUnderCurve = 0.0
	for i in range( 0, numberOfPoints - 1 ):
		Y = psd[ i, 1 ]
		deltaX = np.log10( psd[ i + 1, 0 ] ) - np.log10( psd[ i, 0 ] )
		deltaY = psd[ i + 1, 1 ] - psd[ i, 1 ]
		areaRectangle = Y * deltaX
		areaTriangle = 0.5 * ( deltaX * deltaY )
		areaTotal = areaRectangle + areaTriangle
		areaUnderCurve = areaUnderCurve + areaTotal

	return areaUnderCurve


def getContactNormalsSPAM( labelMap, method='randomWalker', saveData=True, sampleName='', outputDir='', keepPositive='Y'):
	"""Computes the orientations of the inter-particle contacts
	in the scanned data using the spam libraries. This is a loose wrapper.

	The interparticle contacts can be obtained using one of two methods. The
	first method is the itk watershed method. The second method is the random
	walker method. The random walker method is supposed to be more accurate.

	Parameters
	----------
	labelMap : unsigned integer ndarray
		This contains the labelled data

	method : string
		This can be either randomWalker or itkWatershed

	saveData : Bool
		Should the data be saved? Default is True 

	sampleName : string
		Name of the sample. Default is empty (same location as the script)

	outputDir : string
		Where should the data be stored? Default is empty (same location as the script)

	keepPositive : string
		Which axis should be flipped to have it as positive?
		Z - Flips the contact notmals to have all Z as positive
		Y - Flips the contact normals to have all Y as positive
		X - Flips the contact normals to have all X as positive
		Default is Y

	Return
	------
	contTable : float ndarray
		contains details about the interparticle contact orientations.

	"""
	if VERBOSE:
		print( "\nMeasuring contact normals using SPAM library\n" )

	labelledData = Segment.applyPaddingToLabelledMap(labelMap, 2)
	binaryData = np.zeros_like( labelledData )
	binaryData[ np.where( labelledData != 0 ) ] = int( 1 )

	contactVolume, Z, contactsTable, contactingLabels = slab.labelledContacts( labelledData )

	if method == None:
		method = input('Enter the method to use (itk or rw): ')

	if method == 'randomWalker':
		print("\tMeasuring contact using Random Walker\n")
		ortTabSandRW = slab.contacts.contactOrientationsAssembly( labelledData , binaryData , contactingLabels , watershed = "RW" )
		tempOrtsRW = np.zeros_like( ortTabSandRW )
		ortOnlySandRW = ortTabSandRW[ : , 2 : 5 ]

		j = 0
		for i in range( 0 , ortTabSandRW.shape[ 0 ] ):

			if keepPositive == 'Z': axisToKeepPositive = 0
			if keepPositive == 'Y': axisToKeepPositive = 1
			if keepPositive == 'X': axisToKeepPositive = 2

			if ortOnlySandRW[ i , axisToKeepPositive ] < 0:
				ortOnlySandRW[ i ] *= -1

			if ( ortOnlySandRW[ i ] ** 2 ).sum() <= 0.999:
				if VERBOSE: print( "Contact deleted - small contact" )

			else:
				tempOrtsRW[ j ] = ortTabSandRW[ i ]
				j = j + 1

		contactTableRW = tempOrtsRW[ 0 : j , : ]

	elif method == 'itkWatershed':
		if VERBOSE: print( "\tMeasuring contact using ITK watershed\n" )

		ortTabSandITK = slab.contacts.contactOrientationsAssembly( labelledData ,
																	binaryData ,
																	contactingLabels ,
																	watershed="ITK" )

		tempOrtsITK = np.zeros_like( ortTabSandITK )
		ortOnlySandITK = ortTabSandITK[ : , 2 : 5 ]

		j = 0
		for i in range( 0 , ortTabSandITK.shape[ 0 ] ):
			if ortOnlySandITK[ i , 0 ] < 0 : ortOnlySandITK[ i ] *= -1
			if ( ortOnlySandITK[ i ] ** 2 ).sum() <= 0.999 :
				if VERBOSE: print( "Contact deleted - small contact" )

			else:
				tempOrtsITK[ j ] = ortTabSandITK[ i ]
				j = j + 1

		contactTableITK = tempOrtsITK[ 0 : j , : ]

	if method == 'randomWalker' : contTable = contactTableRW
	elif method == 'itkWatershed' : contTable = contactTableITK

	return contTable


def getContactNormals( labelMap, saveData=True, sampleName='', outputDir='' ):
	"""This is a simple implementation that computes the contact normals using
	the random walker algorigthm.

	Parameters
	----------
	labelMap : unsigned int ndarray

	saveData : bool

	sampleName : string

	outputDir : string

	Return
	------
	contactTable : float ndarray
		This table has rows equal to the number of contacts in the sample
		Each contact has columns that contain [0] Particle 1, [1] Particle 2,
		[3] z component of contact normal, [4] y component of contact normal,
		[5] x component of contact normal
	"""


def fabricVariablesSpam( contactTable ):
	"""Computes the fabric tensors using the SPAM library

	Parameter
	---------

	Return
	------

	"""
	orts = contactTable[ :, 2:5]
	F1, F2, F3 = slab.fabricTensor( orts )
	return F1, F2, F3


def fabricVariablesWithUncertainity( contactTable, vectUncert = 0 ):
	"""
	"""
	vectors = contactTable[ :, 2:5]                         # contactTable col 0 is first particle label and col 1 is second particle number
	uncertVectors = vectUncert*(np.ones_like(vectors))      # uncertainity in the vector components

	uncertVectorArray = unp.uarray(vectors,uncertVectors)   # the uncertainity libraries takes vectors and the undertainity of those vectors to make an array

	N = np.zeros((3,3))
	F = np.zeros((3,3))
	Fq = np.zeros((1,1))

	uN = unp.uarray(N,N)
	uF = unp.uarray(F,F)
	uFq = unp.uarray(Fq,Fq)

	for i in range(0,uncertVectorArray.shape[0]):
		uN[0,0] = uN[0,0] + (uncertVectorArray[i,0])*(uncertVectorArray[i,0])
		uN[0,1] = uN[0,1] + (uncertVectorArray[i,0])*(uncertVectorArray[i,1])
		uN[0,2] = uN[0,2] + (uncertVectorArray[i,0])*(uncertVectorArray[i,2])
		uN[1,0] = uN[1,0] + (uncertVectorArray[i,1])*(uncertVectorArray[i,0])
		uN[1,1] = uN[1,1] + (uncertVectorArray[i,1])*(uncertVectorArray[i,1])
		uN[1,2] = uN[1,2] + (uncertVectorArray[i,1])*(uncertVectorArray[i,2])
		uN[2,0] = uN[2,0] + (uncertVectorArray[i,2])*(uncertVectorArray[i,0])
		uN[2,1] = uN[2,1] + (uncertVectorArray[i,2])*(uncertVectorArray[i,1])
		uN[2,2] = uN[2,2] + (uncertVectorArray[i,2])*(uncertVectorArray[i,2])

	uN = uN / uncertVectorArray.shape[0]

	utraceF = uN[0,0] + uN[1,1] + uN[2,2]

	uF[0,0] = 15/2 * ( uN[0,0] - (1/3)*utraceF )
	uF[0,1] = 15/2 * ( uN[0,1] )
	uF[0,2] = 15/2 * ( uN[0,2] )
	uF[1,0] = 15/2 * ( uN[1,0] )
	uF[1,1] = 15/2 * ( uN[1,1] - (1/3)*utraceF )
	uF[1,2] = 15/2 * ( uN[1,2] )
	uF[2,0] = 15/2 * ( uN[2,0] )
	uF[2,1] = 15/2 * ( uN[2,1] )
	uF[2,2] = 15/2 * ( uN[2,2] - (1/3)*utraceF )

	uFq[0,0] = ((3/2)*( uF[0,0]*uF[0,0] + uF[0,1]*uF[0,1] + uF[0,2]*uF[0,2] + uF[1,0]*uF[1,0] + uF[1,1]*uF[1,1] + uF[1,2]*uF[1,2] + uF[2,0]*uF[2,0] + uF[2,1]*uF[2,1] + uF[2,2]*uF[2,2])) ** 0.5

	return uN, uF, uFq


def getCoordinationNumberList( labelledMap, excludeEdgeLabels=True ):
	"""
	"""
	numberOfLabels = labelledMap.max()

	appliedPadding = 2

	coordinationNumberArray = np.zeros( ( numberOfLabels, appliedPadding ) )

	labelledMap = Segment.applyPaddingToLabelledMap(labelledMap, 2)

	for currentLabel in range(1, numberOfLabels + 1):
		print('\nChecking for label ' + str(np.round(currentLabel)))

		contactLabels = slab.contactingLabels( labelledMap, currentLabel, areas=False)
		numberOfContacts = len(contactLabels)

		coordinationNumberArray[currentLabel-1,0] = currentLabel

		edgeLabel = Segment.checkIfEdgeLabel(labelledMap,currentLabel, pad = appliedPadding)

		if edgeLabel == False: coordinationNumberArray[currentLabel-1,1] = numberOfContacts
		if edgeLabel == True: coordinationNumberArray[currentLabel-1,1] = -1

	# labelledMap = Segment.removePaddingFromLabelledMap(padLabMap, 2)

	return coordinationNumberArray


def _getSelectiveCoordinationNumberList(labelMap,labelBoundary=1,outputDir='',sampleName='',saveData=True):
	"""Get coordination numbers of certain labels considering contacts with fixed label ranges
	
	Parameters
	----------
		labelMap: ndarray
		labelBoundary: int
		outputDir: string
		sampleName: string
		saveData: boolean
	
	
	Return
	------
		contactArray: ndArray
	"""
	print('Getting selective coordination number')
	contactArray = np.zeros( ( labelBoundary-1 , 2 ) )
	label = 1

	# For label smaller than labelBoundary
	while label < labelBoundary:
		print('\tChecking label ' + str(np.round(label)) + '/' + str(np.round(labelBoundary-1)))

		# Remove labels that are not the label and not larger than label boundary
		cleanedLabelMap = labelMap
		cleanedLabelMap[ np.where( np.logical_and( (labelMap != label) , (labelMap < labelBoundary) ) ) ] = 0
		contLabels = slab.contacts.contactingLabels( lab=cleanedLabelMap, labels=label, 
													 areas=False, boundingBoxes=None, centresOfMass=None)
		
		numberOfContacts = len(contLabels)
		
		print('\tContacts: ', numberOfContacts)
		
		contactArray[ label-1, 0 ] = label
		contactArray[ label-1, 1 ] = numberOfContacts

		# updated label Map
		label = label + 1

	if saveData == True:
		saveFileName = outputDir + sampleName + '-contactingLabels'
		np.savetxt(saveFileName,contactArray) 

	return contactArray


def computeSphericities(labMap, sampleName='', saveData=True, fixMissingLables=True, outputDir=''):
	"""This function computes sphericities of all the particles in the volume
	
	Sphericity S is defined as:
		S = pi^(1/3) * (6 * Vp)^(2/3) / Ap
	where Vp is the volume of the particle and Ap is the surface area of the particle.
	
	The Volume of the particle is the volume of the all the voxels that make up the particle - thats easy call computeVolumeOfLabel
	
	The surface area of the particle is difficult. Just counting the outermost voxels is no good - can lead to undercounting of the surface areas. 
	The outer surface of the particle can be meshed, using tiny little obedient marching-cubes (https://scikit-image.org/docs/dev/auto_examples/edges/plot_marching_cubes.html).
	The area of the meshed surface can then be computed as the area of mesh elements.
		
	Parameters
	----------
	labMap
	
	sampleName
	
	saveData=True
	
	fixMissingLables
	
	outputDir
	
	Returns
	-------
	sphericityArray: ndArray containing the label number of the particle, the surface area of the particle, the volume of the particle, and the computed sphericity
	
	"""
	if fixMissingLables == True: correctedVolumeToAnalyze = Segment.fixMissingLabels(labMap=labMap, sampleName=sampleName, saveImg=saveData, outputDir=outputDir)
	else: correctedVolumeToAnalyze = labMap
	
	numberOfParticles = correctedVolumeToAnalyze.max()
	sphericityArray = np.zeros( ( numberOfParticles,4 ) )
	
	for ptclNo in range( 1, numberOfParticles + 1 ):
		if VERBOSE: print( 'Checking particle ' + str (ptclNo) + '/' + str(numberOfParticles) )
		areaP = computeSurfaceAreaOfLabel( correctedVolumeToAnalyze, ptclNo )
		volP = computeVolumeOfLabel( correctedVolumeToAnalyze, ptclNo )
		sphericityP = ( ( math.pi )**( 1/3 ) ) * ( ( 6 * volP )**( 2/3 ) ) / areaP
		if VERBOSE: print( '\tA=' + str( np.round( areaP ) ) + ', V=' + str( np.round( volP ) ) +  ', S=' + str( np.round( sphericityP , 2) ) )
		sphericityArray[ ptclNo - 1, 0 ] = ptclNo 	
		sphericityArray[ ptclNo - 1, 1 ] = areaP
		sphericityArray[ ptclNo - 1, 2 ] = volP
		sphericityArray[ ptclNo - 1, 3 ] = sphericityP
		
	if saveData==True: np.savetxt( outputDir + sampleName + '-particleSphericities.csv',sphericityArray, delimiter=',')
	
	return sphericityArray


def computeSurfaceAreaOfLabel(labMap, label):
	"""This function computes the surface area of the particle. It first computes the surface mesh (triagular elements) of 
	the particle with the desired label, then it adds the area of the triangles of the mesh to get the surface area 
	(if there are internal cavities, it will add the areas of those as well)

	Parameters
	----------
	labMap: ndArray of the volume
	label: int label number that needs to be analyzed

	Returns
	-------
	surfaceArea: float with the area of the trangular faces making up the mesh
	"""
	maskedLabel = np.zeros_like(labMap)
	maskedLabel[np.where(labMap == label)] = 1
	
	verts, faces, normals, values = marching_cubes( maskedLabel )
	surfaceArea = mesh_surface_area( verts, faces )
	
	return surfaceArea


def calculateInternalPorosity():
	"""
	"""


def _calculateInternalPorosity(labelMap, label):
	"""
	"""
	
	particle = np.zeros_like(labelMap)
	particle[np.where(labelMap == label)] = 1
	filledParticle=np.zeros_like(labelMap)
	filledParticle[np.where(labelMap == label)] = 1
	
	filledParticle = binary_fill_holes(filledParticle).astype(int)
		
	voidVolume = filledParticle.sum() - particle.sum()
	
	return voidVolume


def calculateEllipsoidOverUnder(labelMap, label, scale=True):
	"""This function measures the relative convexity and concavity of a particle by fitting an ellipsoid to is centroidal axes and
	measuring the volume differences between the two. It produces three values, corresponding to the volume of the particle outside
	the ellipsoid, the amount of empty space inside the ellipsoid, and the overlapped volume between the two shapes.
	
	Parameters
	----------
	labMap: ndArray of the volume
	label: int label number that needs to be analyzed
	scale: bool determining if the volume of the fit ellipsoid will be scaled to match the particle
	
	Returns
	---------
	croppedParticle: ndarray of the particle of interest
	convex: ndarray of the particle volume outside its fit ellipsoid
	concave: ndarray of the empty volume inside the fit ellipsoid
	overlap: ndarray of the particle inside its fit ellipsoid
	convexVol: float with the volume of the particle outside the ellipsoid
	cocaveVol: float with the volume of the empty space inside the ellipsoid
	overlapVol: float with the volume of the particle inside the ellipsoid
	"""
	
	#Repeat of particle normalization from getPrincipalAxesLengths
	
	zyxofLabel = getZYXLocationOfLabel(labelMap,label)
	
	covarianceMatrix = np.cov( zyxofLabel.T )
	
	eigval, eigvec = np.linalg.eig( covarianceMatrix )
	
	meanZ, meanY, meanX = getCenterOfGravityFromZYXLocations( zyxLocationData=zyxofLabel )

	meanMatrix = np.zeros_like( zyxofLabel )

	meanMatrix[ :, 0 ] = meanZ
	meanMatrix[ :, 1 ] = meanY
	meanMatrix[ :, 2 ] = meanX
	
	centeredLocationData = zyxofLabel - meanMatrix
	
	rotationMatrix = np.zeros( ( 3, 3 ) )
	
	rotationMatrix[ :, 0 ] = eigvec[ 0 ]
	rotationMatrix[ :, 1 ] = eigvec[ 1 ]
	rotationMatrix[ :, 2 ] = eigvec[ 2 ]
	
	rotCentPointCloud = ( np.matmul( rotationMatrix, centeredLocationData.T ) ).T #Particle now has 0,0,0 as its mean value and is oriented along orthogonal axes
	
	#Dimensions of particle cloud
	cloudZ = int( rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 0 ].min() )
	cloudY = int( rotCentPointCloud[ :, 1 ].max() - rotCentPointCloud[ :, 1 ].min() )
	cloudX = int( rotCentPointCloud[ :, 2 ].max() - rotCentPointCloud[ :, 2 ].min() )
		
	#Padding chosen to match previous values; arbitrary
	pad = 10
	
	#Create a "big enough" array to fit the particle, with an origin at its center
	croppedParticle = np.zeros((2*(cloudZ+pad), 2*(cloudY+pad), 2*(cloudX+pad)))
	originZ = int(((croppedParticle.shape[0])-1)/2)
	originY = int(((croppedParticle.shape[1])-1)/2)
	originX = int(((croppedParticle.shape[2])-1)/2)
		
	#Create the particle with its center of mass at the origin of the cropped array
	
	for n in np.arange(0,rotCentPointCloud.shape[0]-1):
		croppedParticle[(math.floor( rotCentPointCloud[n,0] )+originZ,math.floor( rotCentPointCloud[n,1] )+originY,math.floor( rotCentPointCloud[n,2] )+originX)]=1
		croppedParticle[(math.ceil( rotCentPointCloud[n,0] )+originZ,math.ceil( rotCentPointCloud[n,1] )+originY,math.ceil( rotCentPointCloud[n,2] )+originX)]=1
	
	croppedParticle = binary_fill_holes(croppedParticle).astype(int)
	
	#Determine axes of ellipsoid based on z,y,z particle axes
	#Optionally, scale axes so particle and ellipsoid volumes are equal
	if scale == True:
		ellVol = math.pi * (4/3) * cloudZ/2 * cloudY/2 * cloudX/2
		partVol = croppedParticle.sum()
		R = (partVol / ellVol)**(1/3)
		a = (cloudZ/2)*R
		b = (cloudY/2)*R
		c = (cloudX/2)*R
	
	elif scale ==False:
		a = (cloudZ/2)
		b = (cloudY/2)
		c = (cloudX/2)
		
	#Generate ellipse with a,b,c dimensions defined by z,y,x axes
	ellipsoid = createVoxelizedEllipsoid( particle = croppedParticle, a = a , b = b, c = c ) 
	
	#Compare ellipse and particle, producing volumes for each measured volume based on value difference between particle and ellipsoid
	overlappedParticle = croppedParticle - ellipsoid
	convex = np.zeros_like(croppedParticle)
	convex[np.where(overlappedParticle == 1)] = 1
	convexVol = convex.sum()
	concave = np.zeros_like(croppedParticle)
	concave[np.where(overlappedParticle == -3)] = 1
	concaveVol = concave.sum()
	overlap = np.zeros_like(croppedParticle)
	overlap[np.where(overlappedParticle == -2)] = 1
	overlapVol = overlap.sum()
	
	return croppedParticle.astype('uint16'), convex.astype('uint16'), concave.astype('uint16'), overlap.astype('uint16'), convexVol, concaveVol, overlapVol


def createVoxelizedEllipsoid( particle, a, b, c ):
	"""This function generates an ellipsoid with its center of mass at the center of particle array and with axes
	a, b, and c corresponding to the z, y, and x centroidal axes of the particle. Axes are oriented along orthogonal coordinates.
	The values of the ellipse voxels are 1.
	
	Parameters
	----------
	particle: ndarray of volume of particle of interest. 
	a: float with the length of the Z axis
	b: float with the length of the Y axis
	c: float with the length of the X axis
	
	Returns
	---------
	ellipsoid: ndarray of volume of fit ellipsoid
	"""
	
	#Create space and origin equivalent to particle of interest
	ellipsoid = np.zeros_like(particle)
	originZ = int(((particle.shape[0])-1)/2)
	originY = int(((particle.shape[1])-1)/2)
	originX = int(((particle.shape[2])-1)/2)
	
	#For each X value along the length of the centroidal axis, create a range of Y values and draw Z as an arc over this range
	xRange = range(int(-c),int(c+1))
	for X in xRange:
		Y = ((1-(X**2)/(c**2))*(b**2))**.5
		yRange = range(int(-Y),int(Y+1))
		for Y in yRange:
			Z = ((1-(int(X)**2)/(c**2)-(int(Y)**2)/(b**2))*(a**2))**.5
			if type(Z) == complex: Z = 0
			if isnan(Z): Z = 0
			Zhi = originZ+int(Z)+1
			Zlo = originZ-int(Z)
			ellipsoid[Zlo:Zhi,int(originY+Y),int(originX+X)]=3
		
	return ellipsoid


def getLargestInscribedSphere(labelMap, label):
	"""This function uses the euclidian distance transform to calculate the distance to the nearest fluid(0) voxel for each voxel in a particle.
	The maximum edt value is then assumed as the radius of the largest possible inscribed sphere for the particle.
	
	Note: Computationally intensive, recommend running cropLabelMap
	
	Parameters
	----------
	labMap: ndarray of the volume
	label: integer of the label of interest
	
	Returns
	inscribedDiameter: scalar of the diameter of the largest sphere which can be inscribed inside the particle
	
	"""
	
	distances=edt(labelMap)
	
	largestSphereRadius = distances.max()
	inscribedDiameter = largestSphereRadius*2
	
	return inscribedDiameter


def getSmallestCircumscribedSphere(labelMap, label):
	"""This function uses zyx coordinates of a particle's surface voxels to calculate the largest distance between voxels in a particle, 
	assuming that this distance is the diameter of the smallest sphere which can circumscribed the particle fully.
	
	Parameters
	----------
	labMap: ndarray of the volume
	label: integer of the label of interest
	
	Returns
	----------
	circumscribedDiameter: scalar of the diameter of the smallest sphere which can be circumscribed around the particle
	
	"""
	#Construct an array of surface voxels, centered at the origin
	
	surfaceVoxels = getSurfaceVoxelsofParticle( labelMap, label, returnWithLabel = True)
	
	zyxofLabel = getZYXLocationOfLabel(surfaceVoxels,label)
	
	meanZ, meanY, meanX = getCenterOfGravityFromZYXLocations( zyxLocationData=zyxofLabel )

	meanMatrix = np.zeros_like( zyxofLabel )

	meanMatrix[ :, 0 ] = meanZ
	meanMatrix[ :, 1 ] = meanY
	meanMatrix[ :, 2 ] = meanX
	
	centeredLocationData = zyxofLabel - meanMatrix
	
	distanceData = np.zeros((centeredLocationData.shape[0],1))
	distanceData[:,0] = (centeredLocationData[:,0]**2 + centeredLocationData[:,1]**2 + centeredLocationData[:,2]**2)**.5
	
	furthestVoxel = int( np.where( (distanceData[:,0] - distanceData[:,0].max()) == 0 )[0].min() )
		
	z = centeredLocationData[furthestVoxel,0]
	y = centeredLocationData[furthestVoxel,1]
	x = centeredLocationData[furthestVoxel,2]
	
	furthestVoxelLocation = np.zeros_like( centeredLocationData )
	furthestVoxelLocation[:,0] = z
	furthestVoxelLocation[:,1] = y
	furthestVoxelLocation[:,2] = x
	
	comparisonMatrix = centeredLocationData - furthestVoxelLocation
	diameters = np.zeros((centeredLocationData.shape[0],1))
		
	diameters = (comparisonMatrix[:,0]**2 + comparisonMatrix[:,1]**2 + comparisonMatrix[:,2]**2)**.5
	
	smallestSphereDiam = diameters.max()
	
	return smallestSphereDiam


def getIrregularityParameter(labelMap, label):
	"""This function the maximum inscribed and minimum circumscribed spheres to calculate a 3D irregularity parameter.
	
	Parameters
	----------
	labMap: ndarray of the volume
	label: integer of the label of interest
	
	Returns
	----------
	irregularityParameter: scalar of the irregularity parameter
	
	"""
	
	inscribed = getLargestInscribedSphere(labelMap, label)
	circumscribed = getSmallestCircumscribedSphere(labelMap, label)
	irregularityParameter = circumscribed / inscribed
	
	return irregularityParameter


def getConvexStatistics(labelMap, label):
	"""This function uses SciPy's convex hull tool to construct a convex hull of a particle and returns its area and volume.
	ConvexHull uses the Qhull library. Info on the computational geometry algorithm is available at www.qhull.org
	
	Parameters
	----------
	labMap: ndarray of the volume
	label: integer of the label of interest
	
	Returns
	----------
	convexArea: scalar of the convex hull's surface area
	convexVolume: scalar of the convex hull's volume
	
	"""
	
	coordinates = getZYXLocationOfLabel(labelMap, label)
	hull = ConvexHull(coordinates)
	convexArea = hull.area
	convexVolume = hull.volume
	
	return convexArea, convexVolume


def getCroppedLabelMap(labelMap, label, pad=0):
	"""This function reduces the label map to the label of interest
	
	Parameters
	----------
	labMap: ndarray of the volume
	label: integer of the label of interest
	pad: bool for optional padding using Segment.applyPaddingToLabelledMap
	
	Returns
	----------
	croppedLabMap: ndarray of the volume of the label map which contains the label of interest and optional padding
		
	"""
	particle = np.zeros_like(labelMap)
	particle[np.where(labelMap == label)]=label
	loc = np.where(particle == label)
	croppedLabMap = particle[ loc[0].min() : loc[0].max()+1,\
							  loc[1].min() : loc[1].max()+1,\
							  loc[2].min() : loc[2].max()+1 ] 
	
	if pad != 0:
		croppedLabMap = Segment.applyPaddingToLabelledMap(croppedLabMap,pad) 
	
	return croppedLabMap


def computeVolumeFractionsOfLabels(labelMap, labelList, saveData=True):
	"""
	"""

#Todo: Test cleaned version
def computeGrainBoundayDistanceStats(referenceMap, gbMap, erodeReference=False, erosionSteps=0,
										saveData=True, dataName='', outputDir='' ):
	"""This code assess how close the boundaries of a grain are to a reference grain

	Parameters
	----------
	referenceMap: The grain boundary reference map
	gbMap: The grain map that is being assessed

	Returns
	-------
	multiMap: The multiplication of the EDM of reference map and gbMap
	flatMultiMap: The 

	"""
	if erodeReference == True:
		referenceMap = erode(referenceMap,iterations=erosionSteps)

	# Normalize values
	invertedReferenceMap = referenceMap // referenceMap.max()
	referenceMap = np.zeros_like(invertedReferenceMap)
	referenceMap[np.where(invertedReferenceMap == 0)] = 1

	gbMap = (gbMap // gbMap.max()).astype('float')
	
	# Conversion done to prevent "too many zeros"
	gbMap[np.where(gbMap == 0)] = np.nan


	# Compute ED map of the reference map
	edmReferenceMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=referenceMap, 
														scaleUp = int(1), 
														saveImg=saveData, 
														sampleName=dataName, 
														outputDir=outputDir)

	# Get multiplication of gb maps with ref maps
	multiMap = (edmReferenceMap * gbMap).astype('float')

	# Get distances
	flatMultiMap = np.ndarray.flatten(multiMap)
	flatMultiMap = flatMultiMap[~np.isnan(flatMultiMap)]

	# Save multiMap and the flattenned map
	if saveData == True:
		np.savetxt( outputDir + dataName + '-GB_DIST.csv', flatMultiMap, delimiter=',')
		tf.imwrite(outputDir + dataName + 'multiMap.tif', multiMap.astype('float'))

	# Compute dist statistics
	minDist = np.nanmin(flatMultiMap)
	maxDist = np.nanmax(flatMultiMap)
	avgDist = np.nanmean(flatMultiMap) 
	sdDist = np.nanstd(flatMultiMap)

	return minDist, maxDist, avgDist, sdDist 

	
def compute2DEqcir(binMap):
	"""Computes the equivalent circular diameter of the binary image
	"""
	binMap = binMap // binMap.max()
	area = binMap.sum()
	diameter = (4*area/math.pi)**0.5

	return diameter


def compute2DCentroid(binMap):
	"""Centroid y and x location
	"""
	props = regionprops(label_image=binMap)
	return props[0].centroid[0], props[0].centroid[1]


def compute2DAR(binMap):
	"""Compute the ratio of the major axes to the minor axis
	"""
	props = regionprops(label_image=binMap)
	ar = props[0].axis_major_length/props[0].axis_minor_length
	return ar 


def computeIOU(binMap1, binMap2):
	"""Simple computation for intersection over union of binary maps
	"""
	binMap1=binMap1//binMap1.max()
	binMap2=binMap2//binMap2.max()

	area1 = binMap1.sum()
	area2 = binMap2.sum()

	intersectionMap = (binMap1+binMap2)//2

	intersection = intersectionMap.sum()
	union = area1 + area2 - intersection

	iou = intersection/union

	return iou