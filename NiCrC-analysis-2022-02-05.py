""" 2022-02-05 NiCrC analysis - EG

"""

# PAC Imports
from pac import Read
from pac import Segment
from pac import Measure
from pac import Plot

# Standard imports
import numpy as np		# Called to create array of subvolume starting points
import matplotlib		# Color map - for orientation plots
import tifffile as tf

#-----------------------------------USER Input start---------------------------------------------*
#------------------------------------------------------------------------------------------------*

# Input and output locations: 
ifl = '/home/eg/Documents/CRG/01-Ni-CrC/NotchScanZoomIn-2022-01-25/'
ofl = '/home/eg/Documents/CRG/01-Ni-CrC/NotchScanZoomIn-2022-01-25/output/particleAnalysis-2022-02-05/'



# FileNames and output file prefix
fileName = 'NotchScanZoomIn-2022-01-25-8bit-650px-CrC-33-132.tif'	# Name of binarized tiff file
dataName = 'CrC-Phase-1-85'											# Prefix for output files - name it smartly
dataFile = ifl + fileName											# The file and location that will be read the bin data

# Calibration
calVal = 0.6197														# mm/vox - keep this as 1 to get sizes in voxel units


# Analysis parameters
edmHVal = 1 													# Minimum peak size when locating local maxima in EDM
radiusRatioVal = 0.85											# Ratio of area radius to smaller particle radius

#------------------------------------------------------------------------------------------------*
#------------------------------------USER Input end----------------------------------------------*


#--------------------------------------------CODE------------------------------------------------*
#------------------------------------------------------------------------------------------------*

print('Running analysis')


subVolumeBinMap = tf.imread(dataFile)//255		# 8 bit file

edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=subVolumeBinMap,
											scaleUp = int(1),
											saveImg=True,
											sampleName=dataName,
											outputDir=ofl )

edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
												h=edmHVal,
												sampleName=dataName,
												saveImg=False,
												outputDir=ofl )

labMap = Segment.segmentUsingWatershed( binaryMapToSeg=subVolumeBinMap,
										edmMapForTopo=edmMap,
										edmPeaksForSeed=edmPeaksMap,
										sampleName=dataName,
										saveImg=True,
										outputDir=ofl,
										addWatershedLine=False)

corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
													pad=int(2),
													areaLimit=700,
													considerEdgeLabels=True,
													checkForSmallParticles=False,
													voxelVolumeThreshold=200,
													radiusCheck=True,
													radiusRatioLimit=radiusRatioVal,
													sampleName=dataName,
													saveImg=True,
													outputDir=ofl )

noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
											pad=0,
											sampleName=dataName,
											saveImg=True,
											outputDir=ofl)

Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis = noEdgeCorLabMap,
								calibrationFactor=calVal,
								saveData=True,
								sampleName=dataName,
								outputDir=ofl )

sph = Measure.computeSphericities(labMap=noEdgeCorLabMap, 
									sampleName=dataName, 
									saveData=True, 
									fixMissingLables=False, 
									outputDir=ofl)

ortsTable = Measure.getPrincipalAxesOrtTable( labelMapForParticleOrientation = noEdgeCorLabMap,
												saveData=True,
												sampleName=dataName,
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
							sampleName=dataName,
							figXSize = 12,
							figYSize = 4.8,
							figFontSize = 15,
							labelName = 'Number of particles')



#------------------------------------------------------------------------------------------------*
#--------------------------------------------CODE------------------------------------------------*
