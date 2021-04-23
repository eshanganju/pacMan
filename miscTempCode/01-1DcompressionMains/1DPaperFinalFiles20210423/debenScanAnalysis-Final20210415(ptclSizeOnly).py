'''
The objective of this code is to assemble the codes in pac directory to:
    1. Read the data files for the deben scans and extract volume
    2. Segment the data (binarize and WS)
    3. Analyze CT data for
        3.1 Correct REV size
        3.2 Particle size for REV size
        3.3 Particle morphology
        3.4 Relative breakage
        3.5 Fabric tensor
    4. Distribution of rel. breakage and fabric tensor
'''
# Importing files
from pac import Read                      # Reads files, cuts into smaller pieces
from pac import Filter                      # Filters files using NLM filter
from pac import Segment                     # Binarizatio and WS
from pac import Measure                     # Calculates particle size, mprphology, contact, breakage
from pac import Plot                        # Plotting functions
import time                                 # The fourth dimension
import matplotlib.pyplot as plt             # matplotlib
import skimage.external.tifffile as tf      # scikit-image
import numpy as np                          # numpy
from uncertainties import unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *

'''
Binarization according to density measurement
Segmentation according to ITK and auto correction
Particle size values (4) and appropriate gradation
Relative breakage according to Einav
Contact according to ITK and RW
Plotting orientations in rose and EAP diagrams
'''
totalTimeStart=time.time()

# run = ['OGF-0N']
# run = [ '2QR-0N-Top','2QR-0N-Middle','2QR-0N-Bottom', '2QR-50N-Middle','2QR-100N-Middle','2QR-500N-Middle', '2QR-1500N-Top', '2QR-1500N-Middle', '2QR-1500N-Bottom']
# run = [ 'OGF-0N', 'OGF-100N', 'OGF-500N', 'OGF-1500N']
run = [ 'OTC-0N' ,'OTC-500N', 'OTC-1500N','OTC-MD-0N', 'OTC-MD-50N', 'OTC-MD-500N', 'OTC-MD-1500N']

for i in run:

	# Void ratios taken from /1DAnalysisFinal20210415/VoidRatioCalculationsCTSamples.ods

	############# 2QR-D #############	
	if i == '2QR-0N-Top' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-0N-Top/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location


	    dataName = i
	    measVoidRatio = 0.734                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 576                                                       # Voxel units - vertical center
	    xCenter = 504                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR-0N-Middle' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-0N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = i
	    measVoidRatio = 0.734                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 427                                                       # Voxel units - vertical center
	    xCenter = 504                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR-0N-Bottom' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-0N-Bottom/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location


	    dataName = i
	    measVoidRatio = 0.734                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 274                                                       # Voxel units - vertical center
	    xCenter = 504                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR-50N-Middle' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-50N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-50N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = i
	    measVoidRatio = 0.726                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 433                                                       # Voxel units - vertical center
	    xCenter = 512                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR-100N-Middle' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-100N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-100N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = i
	    measVoidRatio = 0.722                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 432                                                       # Voxel units - vertical center
	    xCenter = 510                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR-500N-Middle' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-500N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = i
	    measVoidRatio = 0.698                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 438                                                       # Voxel units - vertical center
	    xCenter = 511                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR-1500N-Top' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-1500N-Top/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = i
	    measVoidRatio = 0.591                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 576                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm
	    
	if i == '2QR-1500N-Middle' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-1500N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = i
	    measVoidRatio = 0.591                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 450                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR-1500N-Bottom' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-1500N-Bottom/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = i
	    measVoidRatio = 0.591                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 326                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm
	############# 2QR - D #############	 
	
	############# OGF - D #############
	if i == 'OGF-0N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-0N/' # output folder location ofl
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = i
	    measVoidRatio = 0.635                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.62                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 508                                                       # Voxel units - center of slice
	    yCenter = 437                                                       # Voxel units - vertical center
	    xCenter = 490                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OGF-100N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-100N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-100N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = i
	    measVoidRatio = 0.627                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.62                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 508                                                       # Voxel units - center of slice
	    yCenter = 447                                                       # Voxel units - vertical center
	    xCenter = 490                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OGF-500N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = i
	    measVoidRatio = 0.611                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.62                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 508                                                       # Voxel units - center of slice
	    yCenter = 458                                                       # Voxel units - vertical center
	    xCenter = 490                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OGF-1500N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-1500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = i
	    measVoidRatio = 0.565				                                # Void ratio measured from 1D compression experiment
	    d50 = 0.62                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 508                                                       # Voxel units - center of slice
	    yCenter = 437                                                       # Voxel units - vertical center
	    xCenter = 490                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm
	############# OGF - D #############
	
	############# OTC - D #############
	if i == 'OTC-0N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-0N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = i
	    measVoidRatio = 0.534                                           	# Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 457                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     		# Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OTC-500N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = i
	    measVoidRatio = 0.510                                           # Void ratio measured from 1D compression experiment

	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 455                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     		# Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OTC-1500N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-1500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = i
	    measVoidRatio = 0.492                                           # Void ratio measured from 1D compression experiment

	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 462                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     		# Original GSD
	    eLen = 6*d50          # Edge length in mm
	############# OTC - D #############
	
	############# OTC - MD #############
	if i == 'OTC-MD-0N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-MD-0N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = i
	    measVoidRatio = 0.609                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 405                                                       # Voxel units - vertical center
	    xCenter = 499                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD
	    eLen = 7*d50

	if i == 'OTC-MD-50N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-50N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-MD-50N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = i
	    measVoidRatio = 0.60                                      # Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 486                                                       # Voxel units - vertical center
	    xCenter = 505                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD
	    eLen = 7*d50

	if i == 'OTC-MD-500N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-MD-500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = i
	    measVoidRatio = 0.585                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 495                                                       # Voxel units - vertical center
	    xCenter = 503                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD
	    eLen = 7*d50

	if i == 'OTC-MD-1500N' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-MD-1500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = i
	    measVoidRatio = 0.56                                      # Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 497                                                       # Voxel units - vertical center
	    xCenter = 505                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD
	    eLen = 7*d50
	############# OTC - MD #############

	# Reading and cropping the data file
	# subregionGLIMap = Read.extractSubregionFromTiffSequence(inputFolderLocation,
	#                                                         Z=zCenter,
	#                                                         Y=yCenter,
	#                                                         X=xCenter,
	#                                                         lngt=eLen,
	#                                                         calib=cal,
	#                                                         reference='center',
	#                                                         invImg=False,
	#                                                         saveImg=True, 
	#                                                         outputDir=ofl, 
	#                                                         sampleName=dataName)

	# # Binarization using Otsu
	# binMap = Segment.binarizeAccordingToDensity( gliMapToBinarize=subregionGLIMap,
	#                                              measuredVoidRatio=measVoidRatio,
	#                                              returnThresholdVal=False,
	#                                              saveImg=True,
	#                                              saveData=True,
	#                                              sampleName=dataName,
	#                                              outputDir=ofl)

	# # EDM and particle centers
	# edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
	#                                           scaleUp = int(1),
	#                                           saveImg=False,
	#                                           sampleName=dataName,
	#                                           outputDir=ofl )

	# # EDM peaks
	# edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
	#                                                 h=1,
	#                                                 sampleName=dataName,
	#                                                 saveImg=False,
	#                                                 outputDir=ofl )

	# # Watershed segmentation
	# labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
	#                                         edmMapForTopo=edmMap,
	#                                         edmPeaksForSeed=edmPeaksMap,
	#                                         sampleName=dataName,
	#                                         saveImg=True,
	#                                         outputDir=ofl )

	# # Correction of segmentation errors
	# corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
	#                                              pad=int(2),
	#                                              areaLimit=700,
	#                                              considerEdgeLabels=True,
	#                                              checkForSmallParticles=True,
	#                                              radiusCheck=True,
	#                                              radiusRatioLimit=0.8,
	#                                              sampleName=dataName,
	#                                              saveImg=True,
	#                                              outputDir=ofl )

	# # Removal of edge labels
	# noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
	#                                             pad=0,
	#                                             sampleName=dataName,
	#                                             saveImg=True,
	#                                             outputDir=ofl )
	print('Checking Ferets')
	import tifffile

	noEdgeCorLabMap = tifffile.imread(ofl+dataName+'-noEdgeCorrectedLabelMap.tif')

	# Particle size list
	pss = Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis=noEdgeCorLabMap,
	                                    calibrationFactor=cal,
	                                    saveData=True,
	                                    sampleName=dataName+'(SPAMFeret)',
	                                    outputDir=ofl )

	# Particle size distribution
	psdFeret = Measure.getParticleSizeDistribution( psSummary=pss,
	                                                sizeParam='feretMin',
	                                                saveData=True,
	                                                sampleName=dataName + '(SPAMFeret)',
	                                                outputDir=ofl )

totalTimeEnd = time.time()
totalTimeTaken = totalTimeEnd - totalTimeStart

print('\n\n--------------------------------------**')
print('Total time taken to analyze(mins): ~' + str(totalTimeTaken//60))
