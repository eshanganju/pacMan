""" 2021-09-03 @EG
Main code for the individual segmentation of particles - FGM datasets

This code uses the Segment module to convert the binary data to labelled map

The input needs to be in tiff format, preferably 8 bit
The output will be in tiff as well - labelled map

The particles will be analyzed

"""

#Imports
import tifffile as tf
from pac import Segment
from pac import Measure

# Input and output locations (don't change):
ifl = '/home/chawlahpc2adm/pacInput/particleSegmentationForHamid/'
ofl = '/home/chawlahpc2adm/pacOutput/particleSegmentationForHamid/FGM-final/'

# Change name of the file here
dataFile = ifl + 'segmented2_TIFF.tif'
dataName = 'FGM_CROP_800-1000z_200x200xy_20h_rr08' # Change this to have different output file names

#Reading binary data
binMap = tf.imread(dataFile)[800:1000]

# Getting EDM
edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
	                                   scaleUp = int(1),
	                                   saveImg=True,
	                                   sampleName=dataName,
	                                   outputDir=ofl )

# Getting peaks of EDM	                                        
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
	                                      
noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap, 
					       pad=0, 
					       sampleName=dataName, 
					       saveImg=True, 
					       outputDir=ofl)
	                                
# Particle size anaysis on corlabMap
Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis = noEdgeCorLabMap, 
				calibrationFactor=1, 
				saveData=True, 
				sampleName=dataName, 
				outputDir=ofl)

# Orientation of major axis of particle
Measure.getPrincipalAxesOrtTable( labelMapForParticleOrientation = noEdgeCorLabMap,
				    saveData=True, 
				    sampleName=dataName, 
				    outputDir=ofl)

	                                     

