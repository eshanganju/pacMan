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
import numpy as np

#-----------------------------------USER Input start---------------------------------------------*
#------------------------------------------------------------------------------------------------*

# Input and output locations: 
ifl = '/home/eg/pacInput/fgmInput/'
ofl = '/home/eg/pacOutput/fgmOutput/'

# FileNames and output file prefix
fileName = 'segmented2_cropped.tif'		# Name of binarized tiff file
dataName = 'FGM3-DES_20h_rr08'			# Prefix for output files - name it smartly
dataFile = ifl + fileName				# The file and location that will be read the bin data 

# Calibration
calVal = 0.01135						# mm/vox - keep this as 1 to get sizes in voxel units 

# Subvolume extraction from scan
globalZStart = 100 						# The upper end of the first slice in Z direction 
globalZEnd = 1000						# The lower end of the last slice in Z direction 

csStartX = 100							# Top-left corner X value - check from imageJ
csStartY = 100 							# Top left corner Y value - check from imageJ

subVolumeEdgeLength = 3					# Length (mm units) of the cubical subvolume 

numberOfSubVolumes = 5 					# Number of subvolumes to analyze between gobal Z limits

# Analysis parameters
edmHVal = 2								# Minimum peak size when locating local maxima in EDM
radiusRatioVal = 0.8					# Ratio of area radius to smaller particle radius

#------------------------------------------------------------------------------------------------*
#------------------------------------USER Input end----------------------------------------------*



#--------------------------------------------CODE------------------------------------------------*

# Identify Z values of each slice
upperZFirstSubVolume = globalZStart
upperZLastSubVolume = globalZEnd - ( subVolumeEdgeLength//calVal ) 
upperZSubVolumeList = np.linspace(upperZFirstSubVolume, upperZLastSubVolume, numberOfSubVolumes) 

# Analyze subvolumes corresponding to each z-slice
# This for loop does it one by one - make parallel
for subVolumeNumber in range( 0, numberOfSubVolumes ):
	print( '\n\n\nSubvolume number ' + str( subVolumeNumber + 1 ) )

	# Extract subvolume from the total volume
	subVolumeBinMap = Read.extractSubregionFromTiffSequence( folderDir=dataFile , 
																Z, 
																Y, 
																X, 
																lngt, 
																calib, 
																reference='topLeft',
																invImg=False, 
																saveImg=False, 
																outputDir='', 
																sampleName='XXX')

	# Segment the subvolume using updated watershed + error correction
	#subVolumeEDM = 
	#subVolumeEDMPeaks = 
	#subvolumeLabMap = 
	#subvolumeCorLabMap = 
	#subvolumeNoEdgeCorLabMap = 	

	# Compute and save particle size distributions
	#particleSizeArray = 
	
	# Compute and save aspect ratio distributions
	#aspectRatioSiCParticles = 

	# Compute and save particle orientations
	#sicParticleOrtsArray = 

	# Plot size distribution, aspect ratio, particle orientation

#--------------------------------------------CODE------------------------------------------------*


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
	                                     

