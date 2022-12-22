"""
"""
from pac import Measure
from pac import Segment
import numpy as np
import tifffile as tf
from skimage.morphology import skeletonize
import skimage.measure
from scipy.spatial import ConvexHull
from scipy import ndimage

VERBOSE = True

ofl="/home/eg/Desktop/EG-WateshedAnalysesAvizo/output4-Dilated/"

inputFile = '/home/eg/Desktop/EG-WateshedAnalysesAvizo/EG-WateshedAnalysesAvizo-ALL.tif'

clm = tf.imread(inputFile)

numPtcl = clm.max()

print('Num particles: ' + str(numPtcl))

a = input('-------------')

particleData = np.zeros((numPtcl,6)) 	#Index, Volume, Surface area, NumBranches, Hull area, Hull Volume

genSTL=False
genSkeleton=False
genEDM=False
genEDMSkeleton=False
genPtclVol=False
genSaPtcl=False
genHullData = True


def classifyNodes(particleMap):
	"""Code used to classify nodes as 
			- end node
			- branch node
	"""
	Print('Under construction')

	return nodeClassifiedMap

def numberOfSegments():
	"""
	"""
	Print('Under construction')


def lengthOfSegment():
	"""
	"""
	Print('Under construction')


def tortuosityOfSegment():
	"""
	"""
	Print('Under construction')


def convexHullDataOfParticle(particleMap, dilateParticle=False):
	"""Get convex hull of the particle
	"""
	particleMap = particleMap//particleMap.max()

	if dilateParticle == True:
		particleMap = ndimage.binary_dilation(particleMap,iterations = 1)

	particleData = np.transpose(np.where(particleMap==1))
	
	particleHull = ConvexHull(particleData)

	volume = particleHull.volume
	area = particleHull.area

	return area, volume


for ptclNo in range(1, numPtcl + 1):
	
	print('Checking particle ' + str(ptclNo) + '/' + str(numPtcl))	

	currFileName = 'MgZnSiScan1-' + str(ptclNo)

	particleData[ptclNo-1,0] = ptclNo

	# Extract particle subvolume
	print('\tCropping')
	ptcl = Segment.cropAndPadParticle(labelledMap=clm,
										label=ptclNo,
										pad=5,
										saveData=True,
										fileName= currFileName,
										outputDir=ofl)

	tf.imwrite( (ofl+currFileName+'.tiff'), ptcl.astype('uint8'))

	# Generate STL
	if genSTL == True:
		print('\tMaking Stl')
		Segment._generateInPlaceStlFile( ptcl, 
											stepSize = 1, 
											saveImg=True, 
											sampleName=currFileName, 
											outputDir=ofl)

	# Skeletonize particle subvolume
	if genSkeleton == True:
		print('\tSkeletonizing')
		ptclSkeleton = skeletonize(ptcl)
		ptclSkeleton = ptclSkeleton//ptclSkeleton.max()
		tf.imwrite( (ofl+currFileName+'-skeleton.tiff'), ptclSkeleton.astype('uint8'))

		# Get number of branches
		print('\tGetting number of branches')
		particleData[ptclNo-1,3] = Segment.getNumberOfBranchesOfSkeleton(ptclSkeleton)
		print('\tPtcl ' + str(ptclNo) + ' has ' + str(particleData[numPtcl-1,3]) + 'branches')

	# Compute EDM of particle subvolume
	if genEDM == True:
		print('\tEDMing')
		edmPtcl = Segment.obtainEuclidDistanceMap( binaryMapForEDM=ptcl, 
													scaleUp = int(1), 
													saveImg=False, 
													sampleName=currFileName, 
													outputDir=ofl )

	# Get product of skeleton EDM
	if genEDMSkeleton == True:
		print('\tGetting EDM on skeleton')
		skeletonEDM = ptclSkeleton*edmPtcl

		# Get list of ED along skeleton
		nonZeroEDMVal = (skeletonEDM[np.nonzero(skeletonEDM)]).flatten()
		np.savetxt(ofl+currFileName+'-edmSkeleton.csv', nonZeroEDMVal, delimiter=',')

	# Getting particle volume
	if genPtclVol == True:
		print('\tGetting particle volume')
		particleData[ptclNo-1,1] = np.sum(ptcl)
		print('\t\tVolumee:', str(particleData[ptclNo-1,1]))

	# Surface area of particles
	if genSaPtcl == True:
		print('\tGetting surface area of particle')
		vertices, faces, _, _ = skimage.measure.marching_cubes(volume=ptcl, 
															 		level=None, 
															 		spacing=(1.0, 1.0, 1.0),
																	gradient_direction='descent', 
																	step_size=1, 
																	allow_degenerate=True, 
																	method='lewiner', 
																	mask=None)
		particleData[ptclNo-1,2]= skimage.measure.mesh_surface_area(vertices, faces)
	
	# Hull Data
	if genHullData == True:
		print('\tGetting hull data')
		hullArea, hullVolume = convexHullDataOfParticle(ptcl, dilateParticle=True)
		particleData[ptclNo-1,4] = hullArea
		particleData[ptclNo-1,5] = hullVolume


np.savetxt('MgZnSiScan1-Data.csv',particleData,delimiter=',')