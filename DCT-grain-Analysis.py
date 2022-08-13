"""This pipeline reads the grains and analyzes them for overlap, 2D morphlogy, and grain boundary correspondence
"""


import numpy as np
import tifffile as tf

from pac import Measure

# Input and output locations
iflG = '/home/eg/Desktop/2022-08-02-GrainBoundaryOverlapAnalysis/GrainBoundaryOverlapAnalysis/Comparison-AlignedAndCropped/IndividualGrains-Clean-Solid/'
iflGB = '/home/eg/Desktop/2022-08-02-GrainBoundaryOverlapAnalysis/GrainBoundaryOverlapAnalysis/Comparison-AlignedAndCropped/IndividualGrains-Clean/'
ofl = '/home/eg/Desktop/2022-08-02-GrainBoundaryOverlapAnalysis/GrainBoundaryOverlapAnalysis/Comparison-AlignedAndCropped/AnalysisOutput/'

# Grain properties
# Index(1); Size (3); AR(3); IOU(2); Centroid (3x2); GBDistStats(2x4)
numberOfGrains = 31

for erodeStep in range(2,3):
	print('Erosion: ' + str(erodeStep) )
	GrainDataArray = np.zeros((numberOfGrains,23))

	# Loop over grain number:
	for index in range(0,31):

		print('\tChecking grain ' + str(index))

		# Read grain boundaries from DCT and SEM data 
		semGrainBoundary = tf.imread(iflGB + str(index) + '-SEM-1.tif') 
		cdctGrainBoundary = tf.imread(iflGB + str(index) + '-CDCT2-1.tif')
		hdctGrainBoundary = tf.imread(iflGB + str(index) + '-HDCT-1.tif')
		
		# Index (1)
		GrainDataArray[ index, 0 ] = index
		
		if index !=0 :
			# Read solid grains from DCTs and SEM
			semGrain = tf.imread(iflG + str(index) + '-SEM-1.tif') 
			cdctGrain = tf.imread(iflG + str(index) + '-CDCT2-1.tif')
			hdctGrain = tf.imread(iflG + str(index) + '-HDCT-1.tif')	
			
			#Size (3)
			#SEM
			GrainDataArray[ index, 1 ] = Measure.compute2DEqcir( binMap=semGrain )
			#CDCT
			GrainDataArray[ index, 2 ] = Measure.compute2DEqcir( binMap=cdctGrain )
			#HDCT
			GrainDataArray[ index, 3 ] = Measure.compute2DEqcir( binMap=hdctGrain )
			
			#ARs (3)
			#SEM 
			GrainDataArray[ index, 4 ] = Measure.compute2DAR( binMap=semGrain )
			#CDCT
			GrainDataArray[ index, 5 ] = Measure.compute2DAR( binMap=cdctGrain )
			#HDCT
			GrainDataArray[ index, 6 ] = Measure.compute2DAR( binMap=hdctGrain )
			
			#IOU (2)
			# CDCT
			GrainDataArray[ index, 7 ] = Measure.computeIOU( binMap1=semGrain, 
																binMap2=cdctGrain )
			# HDCT
			GrainDataArray[ index, 8 ] = Measure.computeIOU( binMap1=semGrain, 
																binMap2=hdctGrain )

			#Centroid (3x2 = 6)
			#Centroid-SEM
			centroidLoc1, centroidLoc2 = Measure.compute2DCentroid( binMap=semGrain )
			GrainDataArray[ index, 9 ] = centroidLoc1
			GrainDataArray[ index, 10 ] = centroidLoc2

			#Centroid-CDCT
			centroidLoc1, centroidLoc2 = Measure.compute2DCentroid( binMap=cdctGrain )
			GrainDataArray[ index, 11 ] = centroidLoc1
			GrainDataArray[ index, 12 ] = centroidLoc2
			
			#Centroid-HDCT
			centroidLoc1, centroidLoc2 = Measure.compute2DCentroid( binMap=hdctGrain )
			GrainDataArray[ index, 13 ] = centroidLoc1
			GrainDataArray[ index, 14 ] = centroidLoc2

		#GBDistStats(2x4 = 8)
		# CDCT
		minDistance, maxDistance, avgDistance, sddDistance = Measure.computeGrainBoundayDistanceStats( referenceMap=semGrainBoundary, 
																										gbMap=cdctGrainBoundary,
																										erodeReference=(erodeStep!=0), 
																										erosionSteps=erodeStep,
																										saveData=True, 
																										dataName=str(index)+'-CDCT-Erode-' + str(erodeStep), 
																										outputDir=ofl)
		GrainDataArray[ index, 15 ] = minDistance
		GrainDataArray[ index, 16 ] = maxDistance
		GrainDataArray[ index, 17 ] = avgDistance
		GrainDataArray[ index, 18 ] = sddDistance

		# HDCT
		minDistance, maxDistance, avgDistance, sddDistance = Measure.computeGrainBoundayDistanceStats( referenceMap=semGrainBoundary, 
																										gbMap=hdctGrainBoundary,
																										erodeReference=(erodeStep!=0), 
																										erosionSteps=erodeStep,
																										saveData=True, 
																										dataName=str(index)+'-HDCT-Erode-' + str(erodeStep), 
																										outputDir=ofl)
		GrainDataArray[ index, 19 ] = minDistance
		GrainDataArray[ index, 20 ] = maxDistance
		GrainDataArray[ index, 21 ] = avgDistance
		GrainDataArray[ index, 22 ] = sddDistance


	# Save the outputfile
	np.savetxt( ofl + 'grainDataArray-Erode-' + str(erodeStep) +  '.csv', GrainDataArray, delimiter=",", 
				header="Index,SEM-Dia,CDCT-Dia,HDCT-Dia,SEM-AR,CDCT-AR,HDCT-AR,CDCT-IOU,HDCT-IOU,SEM-CENT-Y,SEM-CENT-X,CDCT-CENT-Y,CDCT-CENT-X,HDCT-CENT-Y,HDCT-CENT-X,CDCT-GBDist-Min,CDCT-GBDist-Max,CDCT-GBDist-Avg,CDCT-GBDist-SD,HDCT-GBDist-Min,HDCT-GBDist-Max,HDCT-GBDist-Avg,HDCT-GBDist-SD")