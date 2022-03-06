""" 2021-09-03 @EG
Main code for the individual segmentation of particles - Hamid's dataset

This code uses the Segment module to convert the binary data to labelled map

The input needs to be in tiff format, preferably 8 bit
The output will be in tiff as well - labelled map



"""

#Imports
import tifffile as tf
from pac import Segment

# Input and output locations (don't change):
ifl = '/home/chawlahpc2adm/pacInput/particleSegmentationForHamid/'
ofl = '/home/chawlahpc2adm/pacOutput/particleSegmentationForHamid/'

# Change name of the file here
dataFile = ifl + 'segmented2.tif'
dataName = 'data1_800_1000_20h_0.8rr' # Change this to have different output file names

#Reading binary data
binMap = tf.imread(dataFile)[800:1000]

# Getting EDM
edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
	                                   scaleUp = int(1),
	                                   saveImg=True,
	                                   sampleName=dataName,
	                                   outputDir=ofl )

# Getting peaks of edm	                                        
edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
	                                         h=2.0, # Decrease this to capture more ptcls
	                                         sampleName=dataName,
	                                         saveImg=True,
	                                         outputDir=ofl )

# Watershed segmentation (without segmentation line)	                                         
labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
	                                 edmMapForTopo=edmMap,
	                                 edmPeaksForSeed=edmPeaksMap,
	                                 sampleName=dataName,
	                                 saveImg=True,
	                                 outputDir=ofl)
	                                 #addWatershedLine=False)

# Correcting watershed segmentation	                                
corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
	                                      pad=int(2),
	                                      areaLimit=700,
	                                      considerEdgeLabels=True,
	                                      checkForSmallParticles=True,
	                                      radiusCheck=True,
	                                      radiusRatioLimit=0.80, # Lower - more aggressive merging
	                                      sampleName=dataName,
	                                      saveImg=True,
	                                      outputDir=ofl )
	                                     

