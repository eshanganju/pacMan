""" 2023-03-04
Segmentation and analysis of the pores within the Adidas doam samples - 2.

"""

# PAC Imports
from pac import Read
from pac import Segment
from pac import Measure
from pac import Plot

# tiffImports
import tifffile as tiffy

# Standard imports
import numpy as np		# Called to create array of subvolume starting points
import matplotlib		# Color map - for orientation plots

#-----------------------------------USER Input start---------------------------------------------*
#------------------------------------------------------------------------------------------------*

# Input and output locations: 
ifl = '/media/chawlahpc2adm/Eshan/EG-AdidasMicrofoam/Output/'
ofl = ifl

# FileNames and output file prefix
fileName = '2_tiffStack-2DMedianFilter1Px-2DOtsuThreshold-Crop430x430x558.tif'	# Name of binarized tiff file
dataName = '2_Adidas-Microfoam-Hval2-rr08'													# Prefix for output files - name it smartly
dataFile = ifl + fileName														# The file and location that will be read the bin data

# Calibration
calVal = 10.79																	# mm/vox - keep this as 1 to get sizes in voxel units

# Analysis parameters
edmHVal = 2 																	# Minimum peak size when locating local maxima in EDM
radiusRatioVal = 0.8															# Ratio of area radius to smaller particle radius

#------------------------------------------------------------------------------------------------*
#------------------------------------USER Input end----------------------------------------------*


#--------------------------------------------CODE------------------------------------------------*
#------------------------------------------------------------------------------------------------*

print('Running analysis')

# Bimodal distribution allows us to use OTSU's method
binMap = tiffy.imread(dataFile)

edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
											scaleUp = int(1),
											saveImg=True,
											sampleName=dataName,
											outputDir=ofl )

edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
												h=edmHVal,
												sampleName=dataName,
												saveImg=False,
												outputDir=ofl )

labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
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
													voxelVolumeThreshold=10,
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
