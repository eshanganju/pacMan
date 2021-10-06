""" 2021-10-06 SiC paper code
This code takes the Segmented (binarized) data from the DES scan or the ML code and then:
- Extract subvolumes from the binarized dataset
- Segments the dataset to get individual particles
- Computes the particle size table
- Computes the particle orientation
- Plots the particle orientation
- Plot particle size distribution
- Plot particle aspect ratio
"""

#Imports
import tifffile as tf	# read the binary file
from pac import Segment
from pac import Measure
from pac import Plot

# Input and output locations (don't change):
ifl = '/home/chawlahpc2adm/pacInput/particleSegmentationForHamid/'
ofl = '/home/chawlahpc2adm/pacOutput/particleSegmentationForHamid/newAnalysis-2021-10-04/'

# FileNames and output file prefix
dataFile = ifl + 'segmented2_cropped.tif'
dataName = 'FGM3-ML_200x200x200zyx_20h_rr08' 


for subVolumes in volumes:
	# Extract subvolume from the total volume

	# Segment the subvolume using updated watershed + error correction

	# Compute and save particle size distributions

	# Compute and save aspect ratio distributions

	# Compute and save particle orientations

	# Plot size distribution, aspect ratio, particle orientation



##--------------------------Extra code --------------------------------------------#
##Reading binary data
#binMap = tf.imread(dataFile)[800:1000]
#
## Getting EDM
#edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
#	                                   scaleUp = int(1),
#	                                   saveImg=True,
#	                                   sampleName=dataName,
#	                                   outputDir=ofl )
#
## Getting peaks of EDM	                                        
#edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
#	                                         h=2, # Decrease this to capture more ptcls
#	                                         sampleName=dataName,
#	                                         saveImg=True,
#	                                         outputDir=ofl )
#
## Watershed segmentation (without segmentation line)	                                         
#labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
#	                                 edmMapForTopo=edmMap,
#	                                 edmPeaksForSeed=edmPeaksMap,
#	                                 sampleName=dataName,
#	                                 saveImg=True,
#	                                 outputDir=ofl)
#	                                 #addWatershedLine=False)
#
## Correcting watershed segmentation	                                
#corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
#	                                      pad=int(2),
#	                                      areaLimit=700,
#	                                      considerEdgeLabels=True,
#	                                      checkForSmallParticles=True,
#	                                      radiusCheck=True,
#	                                      radiusRatioLimit=0.80, # Lower - more aggressive merging
#	                                      sampleName=dataName,
#	                                      saveImg=True,
#	                                      outputDir=ofl )
#	                                      
#noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap, 
#					       pad=0, 
#					       sampleName=dataName, 
#					       saveImg=True, 
#					       outputDir=ofl)
#	                                
## Particle size anaysis on corlabMap
#Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis = noEdgeCorLabMap, 
#				calibrationFactor=1,	# mm per voxel 
#				saveData=True, 
#				sampleName=dataName, 
#				outputDir=ofl)
#
## Orientation of major axis of particle
#ortsTable = Measure.getPrincipalAxesOrtTable( labelMapForParticleOrientation = noEdgeCorLabMap,
#						saveData=True, 
#					    	sampleName=dataName, 
#						outputDir=ofl)
#
## Plot orientations
#Plot.plotOrientationsSPAM(ortsTable[:,1:],
#                         projection="lambert",
#                         plot="both",
#                         binValueMin=0,
#                         binValueMax=20,
#                         binNormalisation = False,
#                         numberOfRings = 9,
#                         # the size of the dots in the points plot (5 OK for many points, 25 good for few points/debugging)
#                         pointMarkerSize = 8,
#                         cmap = matplotlib.pyplot.cm.Reds,
#                         title = "",
#                         subtitle = {"points":"","bins":""},
#                         saveFigPath = ofl,
#                         sampleName=dataName,
#                         figXSize = 12,
#                         figYSize = 4.8,
#                         figFontSize = 15,
#                         labelName = 'Number of particles')
#
	                                     

