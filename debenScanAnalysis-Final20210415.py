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



# run=['2QR_D_0_Top','2QR_D_0_Middle','2QR_D_0_Bottom','2QR_D_50','2QR_D_100','2QR_D_500', '2QR_D_1500_Top', '2QR_D_1500_Middle', '2QR_D_1500_Bottom']

# run=['OGF_D_0', 'OGF_D_100', 'OGF_D_500', 'OGF_D_1500']

run=['OTC_D_0' ,'OTC_D_500', 'OTC_D_1500', 'OTC_MD_0', 'OTC_MD_50', 'OTC_MD_500', 'OTC_MD_1500' ]

for i in run:

	# Void ratios taken from /1DAnalysisFinal20210415/VoidRatioCalculationsCTSamples.ods

	############# 2QR-D #############	
	if i == '2QR_D_0_Top' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-0N-Top/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location


	    dataName = '2qr-0N-Top'
	    measVoidRatio = 0.734                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 576                                                       # Voxel units - vertical center
	    xCenter = 504                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR_D_0_Middle' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-0N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = '2QR-0N-Middle'
	    measVoidRatio = 0.734                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 427                                                       # Voxel units - vertical center
	    xCenter = 504                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR_D_0_Bottom' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-0N-Bottom/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location


	    dataName = '2qr-0N-Bottom'
	    measVoidRatio = 0.734                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 274                                                       # Voxel units - vertical center
	    xCenter = 504                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR_D_50' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-50N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-50N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = '2qr-50N-Middle'
	    measVoidRatio = 0.726                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011932                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 433                                                       # Voxel units - vertical center
	    xCenter = 512                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR_D_100' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-100N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-100N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = '2qr-100N-Middle'
	    measVoidRatio = 0.722                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 432                                                       # Voxel units - vertical center
	    xCenter = 510                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR_D_500' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-500N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = '2qr-500N-2-Middle'
	    measVoidRatio = 0.698                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 438                                                       # Voxel units - vertical center
	    xCenter = 511                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR_D_1500_Top' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-1500N-Top/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = '2qr-1500N-Top'
	    measVoidRatio = 0.591                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 576                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm
	    
	if i == '2QR_D_1500_Middle' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-1500N-Middle/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = '2qr-1500N-Middle'
	    measVoidRatio = 0.591                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.73                                                          # D50 in mm - original gradation
	    cal = 0.011931                                                      # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 450                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == '2QR_D_1500_Bottom' :
	    inputFolderLocation = '/home/eg/codes/pacInput/2QR-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/2QR-1500N-Bottom/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

	    dataName = '2qr-1500N-Bottom'
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
	if i == 'OGF_D_0' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-0N/' # output folder location ofl
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = 'ogf-0N'
	    measVoidRatio = 0.635                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.62                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 508                                                       # Voxel units - center of slice
	    yCenter = 437                                                       # Voxel units - vertical center
	    xCenter = 490                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OGF_D_100' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-100N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-100N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = 'ogf-100N'
	    measVoidRatio = 0.627                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.62                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 508                                                       # Voxel units - center of slice
	    yCenter = 447                                                       # Voxel units - vertical center
	    xCenter = 490                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OGF_D_500' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = 'ogf-500N'
	    measVoidRatio = 0.611                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.62                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 508                                                       # Voxel units - center of slice
	    yCenter = 458                                                       # Voxel units - vertical center
	    xCenter = 490                                                       # Voxel units - horizontal center
	    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OGF_D_1500' :
	    print('1500!!!!')
	    inputFolderLocation = '/home/eg/codes/pacInput/OGF-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OGF-1500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = 'ogf-1500N'
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
	if i == 'OTC_D_0' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-0N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-0N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = 'otc-0N'
	    measVoidRatio = 0.534                                           	# Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 457                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     		# Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OTC_D_500' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = 'otc-500N'
	    measVoidRatio = 0.510                                           # Void ratio measured from 1D compression experiment

	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 455                                                       # Voxel units - vertical center
	    xCenter = 507                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     		# Original GSD
	    eLen = 6*d50          # Edge length in mm

	if i == 'OTC_D_1500' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-1500N/'
	    ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/OTC-1500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = 'otc-1500N'
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
	if i == 'OTC_MD_0' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-0N/'
	    ofl = '/home/eg/codes/pacOutput/OTC-MD-0N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details 0N:
	    dataName = 'otc-MD-0N'
	    measVoidRatio = 0.609                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 405                                                       # Voxel units - vertical center
	    xCenter = 499                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD
	    eLen = 7*d50

	if i == 'OTC_MD_50' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-50N/'
	    ofl = '/home/eg/codes/pacOutput/OTC-MD-50N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = 'otc-MD-50N'
	    measVoidRatio = 0.60                                      # Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 486                                                       # Voxel units - vertical center
	    xCenter = 505                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD
	    eLen = 7*d50

	if i == 'OTC_MD_500' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-500N/'
	    ofl = '/home/eg/codes/pacOutput/OTC-MD-500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = 'otc-MD-500N'
	    measVoidRatio = 0.585                                     # Void ratio measured from 1D compression experiment
	    d50 = 0.72                                                          # D50 in mm - original gradation
	    cal = 0.01193                                                       # calibration from CT mm/voxel
	    zCenter = 513                                                       # Voxel units - center of slice
	    yCenter = 495                                                       # Voxel units - vertical center
	    xCenter = 503                                                       # Voxel units - horizontal center
	    origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD
	    eLen = 7*d50

	if i == 'OTC_MD_1500' :
	    inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-1500N/'
	    ofl = '/home/eg/codes/pacOutput/OTC-MD-1500N/'
	    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

	    # Data details:
	    dataName = 'otc-MD-1500N'
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
	subregionGLIMap = Read.extractSubregionFromTiffSequence(inputFolderLocation,
	                                                        Z=zCenter,
	                                                        Y=yCenter,
	                                                        X=xCenter,
	                                                        lngt=eLen,
	                                                        calib=cal,
	                                                        reference='center',
	                                                        invImg=False,
	                                                        saveImg=True, 
	                                                        outputDir=ofl, 
	                                                        sampleName=dataName)

	# Binarization using Otsu
	binMap = Segment.binarizeAccordingToDensity( gliMapToBinarize=subregionGLIMap,
	                                             measuredVoidRatio=measVoidRatio,
	                                             returnThresholdVal=False,
	                                             saveImg=True,
	                                             saveData=True,
	                                             sampleName=dataName,
	                                             outputDir=ofl)

	# EDM and particle centers
	edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
	                                          scaleUp = int(1),
	                                          saveImg=False,
	                                          sampleName=dataName,
	                                          outputDir=ofl )

	# EDM peaks
	edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
	                                                h=1,
	                                                sampleName=dataName,
	                                                saveImg=False,
	                                                outputDir=ofl )

	# Watershed segmentation
	labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
	                                        edmMapForTopo=edmMap,
	                                        edmPeaksForSeed=edmPeaksMap,
	                                        sampleName=dataName,
	                                        saveImg=True,
	                                        outputDir=ofl )

	# Correction of segmentation errors
	corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
	                                             pad=int(2),
	                                             areaLimit=700,
	                                             considerEdgeLabels=True,
	                                             checkForSmallParticles=False,
	                                             radiusCheck=True,
	                                             radiusRatioLimit=0.8,
	                                             sampleName=dataName,
	                                             saveImg=True,
	                                             outputDir=ofl )

	# Removal of edge labels
	noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
	                                            pad=0,
	                                            sampleName=dataName,
	                                            saveImg=True,
	                                            outputDir=ofl )

	# Particle size list
	pss = Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis=noEdgeCorLabMap,
	                                    calibrationFactor=cal,
	                                    saveData=True,
	                                    sampleName=dataName,
	                                    outputDir=ofl )

	# Particle size distribution
	psdFeret = Measure.getParticleSizeDistribution( psSummary=pss,
	                                                sizeParam='feretMin',
	                                                saveData=True,
	                                                sampleName=dataName,
	                                                outputDir=ofl )

totalTimeEnd = time.time()
totalTimeTaken = totalTimeEnd - totalTimeStart

print('\n\n--------------------------------------**')
print('Total time taken to analyze(mins): ~' + str(totalTimeTaken//60))
