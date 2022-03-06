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
ifl = '/home/chawlahpc2adm/Desktop/EG-NiCrC/'
ofl = '/home/chawlahpc2adm/Desktop/EG-NiCrC/output/'



# FileNames and output file prefix
fileName = 'NotchScanZoomIn-2022-01-25-8bit-650px-Void-0-32.tif'	# Name of binarized tiff file
dataName = 'Void-Phase-2-60-Center'									# Prefix for output files - name it smartly
dataFile = ifl + fileName											# The file and location that will be read the bin data

# Subvolume extraction from scan
globalZStart = 150												# The upper end of the first slice in Z direction
globalZEnd = 500												# The lower end of the last slice in Z direction

csStartX = 145													# Top-left corner X value - check from imageJ
csStartY = 145 													# Top left corner Y value - check from imageJ

subVolumeEdgeLength = 350										# Length (calibrated units) of the cubical subvolume

# Calibration
calVal = 1														# mm/vox - keep this as 1 to get sizes in voxel units


# Analysis parameters
edmHVal = 2 													# Minimum peak size when locating local maxima in EDM
radiusRatioVal = 0.60											# Ratio of area radius to smaller particle radius

#------------------------------------------------------------------------------------------------*
#------------------------------------USER Input end----------------------------------------------*


#--------------------------------------------CODE------------------------------------------------*
#------------------------------------------------------------------------------------------------*

print('Running analysis')


#subVolumeBinMap = tf.imread(dataFile)//255		# 8 bit file

subVolumeBinMap = Read.extractSubregionFromTiffFile(fileDir=dataFile,
													Z=globalZStart,
													Y=csStartY,
													X=csStartX,
													lngt=subVolumeEdgeLength,
													calib=calVal,
													zReference='low',
													xyReference='topLeft',
													invImg=False,
													saveImg=True,
													outputDir=ofl,
													sampleName=dataName)

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
													checkForSmallParticles=True,
													voxelVolumeThreshold=125,
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
