""" 2021-10-06 SiC paper code - EG
This code takes the Segmented (binarized) data from the DES scan or the ML code and then:
- Extract subvolumes from the binarized dataset
- Segments the dataset to get individual particles
- Computes the particle size table
- Computes the particle orientation
- Plots the particle orientation
- Plot particle size distribution
- Plot particle aspect ratio
"""

# PAC Imports
from pac import Read
from pac import Segment
from pac import Measure
from pac import Plot

# Standard imports
import numpy as np		# Called to create array of subvolume starting points
import matplotlib		# Color map - for orientation plots

#-----------------------------------USER Input start---------------------------------------------*
#------------------------------------------------------------------------------------------------*

# Input and output locations: 
ifl = '/home/chawlahpc2adm/pacInput/particleSegmentationForHamid/'
ofl = '/home/chawlahpc2adm/pacOutput/FGM-3_DES/1_3umDataset/'



# FileNames and output file prefix
fileName = 'segmented2.tif'										# Name of binarized tiff file
dataName = 'FGM3-ML_1_3um_10h_rr08'								# Prefix for output files - name it smartly
dataFile = ifl + fileName										# The file and location that will be read the bin data

# Calibration
calVal = 1														# mm/vox - keep this as 1 to get sizes in voxel units

# Subvolume extraction from scan
globalZStart = 40 												# The upper end of the first slice in Z direction
globalZEnd = 2170												# The lower end of the last slice in Z direction

csStartX = 0													# Top-left corner X value - check from imageJ
csStartY = 105 													# Top left corner Y value - check from imageJ

subVolumeEdgeLength = 288										# Length (calibrated units) of the cubical subvolume

numberOfSubVolumes = 5 											# Number of subvolumes to analyze between gobal Z limits

# Analysis parameters
edmHVal = 1														# Minimum peak size when locating local maxima in EDM
radiusRatioVal = 0.8											# Ratio of area radius to smaller particle radius

#------------------------------------------------------------------------------------------------*
#------------------------------------USER Input end----------------------------------------------*


#--------------------------------------------CODE------------------------------------------------*
#------------------------------------------------------------------------------------------------*

# Identify Z values of each slice
upperZFirstSubVolume = globalZStart
upperZLastSubVolume = globalZEnd - ( subVolumeEdgeLength//calVal ) 
upperZSubVolumeList = np.linspace(upperZFirstSubVolume, upperZLastSubVolume, numberOfSubVolumes) 

# Analyze subvolumes corresponding to each z-slice

#TODO: Parallelize this loop  - check joblib - https://stackoverflow.com/questions/9786102/how-do-i-parallelize-a-simple-python-loop

for subVolumeNumber in range( 1, numberOfSubVolumes ):
	print( '\n\n\nSubvolume number ' + str( subVolumeNumber + 1 ) + '/' + str(numberOfSubVolumes) )

	zSlice = upperZSubVolumeList[subVolumeNumber]

	subVolumeBinMap = Read.extractSubregionFromTiffFile(fileDir=dataFile,
														Z=zSlice,
														Y=csStartY,
														X=csStartX,
														lngt=subVolumeEdgeLength,
														calib=calVal,
														zReference='low',
														xyReference='topLeft',
														invImg=False,
														saveImg=False,
														outputDir=ofl,
														sampleName=(dataName+'-'+str(subVolumeNumber)))

	edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=subVolumeBinMap,
												scaleUp = int(1),
												saveImg=False,
												sampleName=dataName+'-'+str(subVolumeNumber),
												outputDir=ofl )

	edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
													h=edmHVal,
													sampleName=dataName+'-'+str(subVolumeNumber),
													saveImg=False,
													outputDir=ofl )

	labMap = Segment.segmentUsingWatershed( binaryMapToSeg=subVolumeBinMap,
											edmMapForTopo=edmMap,
											edmPeaksForSeed=edmPeaksMap,
											sampleName=dataName+'-'+str(subVolumeNumber),
											saveImg=True,
											outputDir=ofl,
											addWatershedLine=False)

	corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
														pad=int(2),
														areaLimit=700,
														considerEdgeLabels=True,
														checkForSmallParticles=True,
														radiusCheck=True,
														radiusRatioLimit=radiusRatioVal,
														sampleName=dataName+'-'+str(subVolumeNumber),
														saveImg=True,
														outputDir=ofl )

	noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
												pad=0,
												sampleName=dataName+'-'+str(subVolumeNumber),
												saveImg=True,
												outputDir=ofl)

	Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis = noEdgeCorLabMap,
									calibrationFactor=calVal,
									saveData=True,
									sampleName=dataName+'-'+str(subVolumeNumber),
									outputDir=ofl )

	ortsTable = Measure.getPrincipalAxesOrtTable( labelMapForParticleOrientation = noEdgeCorLabMap,
													saveData=True,
													sampleName=dataName+'-'+str(subVolumeNumber),
													outputDir=ofl)

	Plot.plotOrientationsSPAM( ortsTable[:,1:],
								projection="lambert",
								plot="both",
								binValueMin=0,
								binValueMax=20,
								binNormalisation=False,
								numberOfRings=9,
								pointMarkerSize=8,
								cmap=matplotlib.pyplot.cm.Reds,
								title="",
								subtitle={"points":"","bins":""},
								saveFigPath=ofl,
								sampleName=dataName+'-'+str(subVolumeNumber),
								figXSize = 12,
								figYSize = 4.8,
								figFontSize = 15,
								labelName = 'Number of particles')

#------------------------------------------------------------------------------------------------*
#--------------------------------------------CODE------------------------------------------------*
