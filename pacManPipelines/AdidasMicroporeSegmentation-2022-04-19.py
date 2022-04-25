""" 2023-03-04
Segmentation and analysis of the pores within the Adidas doam samples.

"""

# PAC Imports
from pac import Read
from pac import Segment
from pac import Measure
from pac import Plot

# Standard imports
import numpy as np		# Called to create array of subvolume starting points
import matplotlib		# Color map - for orientation plots
import tifffile as tf 	# image import

#-----------------------------------USER Input start---------------------------------------------*
#------------------------------------------------------------------------------------------------*

# Input and output locations: 
ifl = '/home/chawlahpc2adm/Desktop/eshanganju/'
ofl = '/home/chawlahpc2adm/Desktop/eshanganju/outPut/'

# FileNames and output file prefix
fileName = 'microFoam_1_Bin_crop-fillHoles.tif'					# Name of binarized tiff file
dataName = 'microFoam1_H3_RR80'									# Prefix for output files - name it smartly
dataFile = ifl + fileName										# The file and location that will be read the bin data

# Calibration
calVal = 1														# mm/vox - keep this as 1 to get sizes in voxel units

# Analysis parameters
edmHVal = 3 													# Minimum peak size when locating local maxima in EDM
radiusRatioVal = 0.8											# Ratio of area radius to smaller particle radius

#------------------------------------------------------------------------------------------------*
#------------------------------------USER Input end----------------------------------------------*


#--------------------------------------------CODE------------------------------------------------*
#------------------------------------------------------------------------------------------------*

print('Running analysis')


# subVolumeBinMap = Read.extractSubregionFromTiffFile(fileDir=dataFile,
# 													Z=globalZStart,
# 													Y=csStartY,
# 													X=csStartX,
# 													lngt=subVolumeEdgeLength,
# 													calib=calVal,
# 													zReference='low',
# 													xyReference='topLeft',
# 													invImg=False,
# 													saveImg=True,
# 													outputDir=ofl,
# 													sampleName=dataName)

# subVolumeBinMap = tf.imread(dataFile)//255

# edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=subVolumeBinMap,
# 											scaleUp = int(1),
# 											saveImg=True,
# 											sampleName=dataName,
# 											outputDir=ofl )


# edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
# 												h=edmHVal,
# 												sampleName=dataName,
# 												saveImg=True,
# 												outputDir=ofl )

# labMap = Segment.segmentUsingWatershed( binaryMapToSeg=subVolumeBinMap,
# 										edmMapForTopo=edmMap,
# 										edmPeaksForSeed=edmPeaksMap,
# 										sampleName=dataName,
# 										saveImg=True,
# 										outputDir=ofl,
# 										addWatershedLine=False)

# corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
# 													pad=int(2),
# 													areaLimit=700,
# 													considerEdgeLabels=True,
# 													checkForSmallParticles=True,
# 													voxelVolumeThreshold=10,
# 													radiusCheck=True,
# 													radiusRatioLimit=radiusRatioVal,
# 													sampleName=dataName,
# 													saveImg=True,
# 													outputDir=ofl )

# noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
# 											pad=0,
# 											sampleName=dataName,
# 											saveImg=True,
# 											outputDir=ofl)

noEdgeFileNameToRead = ofl + dataName + '-noEdgeCorrectedLabelMap.tif'
noEdgeCorLabMap = tf.imread(noEdgeFileNameToRead)

Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis = noEdgeCorLabMap,
								getCaDia=True,
								getFeretDia=False,
								calibrationFactor=calVal,
								saveData=True,
								sampleName=dataName,
								outputDir=ofl )

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
