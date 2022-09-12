"""
Pipeline to segment pores a and solid phase
"""

from pac import Segment
from pac import Measure
import tifffile as tf
import numpy as np

# Read
binMap = tf.imread("/home/eg/Desktop/Al-AM-samples/output/sample2-Solid-Bin150gli-CROP.tif")
binMap = binMap//binMap.max()

dataName = "Sample-2-Solid-"
lofMap = np.zeros_like(binMap)

# Parameters
# ----------

# Calibration
calVal = 1														# mm/vox - keep this as 1 to get sizes in voxel units

# Analysis parameters
edmHVal = 1 													# Minimum peak size when locating local maxima in EDM
radiusRatioVal = 0.8											# Ratio of area radius to smaller particle radius

# Output locations
ofl = "/home/eg/Desktop/Al-AM-samples/output/pacManAnalysis/"

# SizeLimit
limitMaxInscSphere = 100

# Segmentation
# ------------
edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
											scaleUp = int(1),
											saveImg=True,
											sampleName=dataName,
											outputDir=ofl )

edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
												h=edmHVal,
												sampleName=dataName,
												saveImg=True,
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
													checkForSmallParticles=False,
													voxelVolumeThreshold=10,
													radiusCheck=True,
													radiusRatioLimit=radiusRatioVal,
													sampleName=dataName,
													saveImg=True,
													outputDir=ofl )

# Size analysis
# -------------

# Get edm of binarized map

# Check the max incribed circle raidus of each label
# for label in range(1,corLabMap.max()+1):
	# in EDM map, get max value where label is label in corlabMap

	# if max edm is smaller than threshold, add it to LOF map


# Cleanup - Adjacent occurrence of small pores and particles
# -------
