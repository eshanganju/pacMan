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
from pac import Plot

# Input and output locations (don't change):
ifl = '/home/eg/pacInput/DSCoVERAnalysis_FGM3_2021-10-03/2021-10-03-Analysis-Final/'
ofl = '/home/eg/pacInput/DSCoVERAnalysis_FGM3_2021-10-03/2021-10-03-Analysis-Final/output/'

# Change name of the file here
dataFile = ifl + 'SiC_2_CMASK_crop_final_volfrac_100.tif'
dataName = 'FGM3-DSC0VER_100x100x100zyx_10h_rr08' # Change this to have different output file names

#Reading binary data
binMap = tf.imread(dataFile)	# Makes it 0 and 1

# Getting EDM
edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
	                                   scaleUp = int(1),
	                                   saveImg=True,
	                                   sampleName=dataName,
	                                   outputDir=ofl )

# Getting peaks of EDM	                                        
edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
	                                         h=1.0, # Decrease this to capture more ptcls
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
Measure.getParticleSizeArraye( labelledMapForParticleSizeAnalysis = noEdgeCorLabMap, 
				calibrationFactor=1,	# mm per voxel 
				saveData=True, 
				sampleName=dataName, 
				outputDir=ofl)

# Orientation of major axis of particle
ortsTable = Measure.getPrincipalAxesOrtTable( labelMapForParticleOrientation = noEdgeCorLabMap,
						saveData=True, 
					    	sampleName=dataName, 
						outputDir=ofl)

# Plot orientations
Plot.plotOrientationsSPAM(ortTable[:,1:],
                         projection="lambert",
                         plot="both",
                         binValueMin=0,
                         binValueMax=10,
                         binNormalisation = False,
                         numberOfRings = 9,
                         # the size of the dots in the points plot (5 OK for many points, 25 good for few points/debugging)
                         pointMarkerSize = 8,
                         cmap = matplotlib.pyplot.cm.RdBu_r,
                         title = "",
                         subtitle = {"points":"","bins":""},
                         saveFigPath = ofl,
                         sampleName=dataName,
                         figXSize = 6.1,
                         figYSize = 4.8,
                         figFontSize = 15,
                         labelName = 'Number of particles')

	                                     

