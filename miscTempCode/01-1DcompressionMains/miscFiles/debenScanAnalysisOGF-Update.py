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
from pac import Read                        # Reads files, cuts into smaller pieces
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

# 0 (0 MPa), 100 (2 MPa), 500 (10 MPa), 1500 (30 MPa), 4500 (90 MPa)

for i in [1500]:
    
    print('Running i = ' + str(i))
    
    if i == 0 :
        inputFolderLocation = '/home/eg/codes/pacInput/OGF-0N/'
        ofl = '/home/eg/codes/pacOutput/OGF-0N/' # output folder location ofl
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

    if i == 100 :
        inputFolderLocation = '/home/eg/codes/pacInput/OGF-100N/'
        ofl = '/home/eg/codes/pacOutput/OGF-100N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

        # Data details 0N:
        dataName = 'ogf-100N'
        measVoidRatio = 0.625                                     # Void ratio measured from 1D compression experiment
        d50 = 0.62                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 508                                                       # Voxel units - center of slice
        yCenter = 447                                                       # Voxel units - vertical center
        xCenter = 490                                                       # Voxel units - horizontal center
        origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD

    if i == 500 :
        inputFolderLocation = '/home/eg/codes/pacInput/OGF-500N/'
        ofl = '/home/eg/codes/pacOutput/OGF-500N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

        # Data details 0N:
        dataName = 'ogf-500N'
        measVoidRatio = 0.615                                     # Void ratio measured from 1D compression experiment
        d50 = 0.62                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 508                                                       # Voxel units - center of slice
        yCenter = 458                                                       # Voxel units - vertical center
        xCenter = 490                                                       # Voxel units - horizontal center
        origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD

    if i == 1500 :
        print('1500!!!!')
        inputFolderLocation = '/home/eg/codes/pacInput/OGF-1500N/'
        ofl = '/home/eg/codes/pacOutput/OGF-1500N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

        # Data details 0N:
        dataName = 'ogf-1500N'
        measVoidRatio = 0.571  #0.565                             # Void ratio measured from 1D compression experiment
        d50 = 0.62                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 508                                                       # Voxel units - center of slice
        yCenter = 437                                                       # Voxel units - vertical center
        xCenter = 490                                                       # Voxel units - horizontal center
        origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD

    eLen = 6*d50          # Edge length in mm

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
                                                    h=3,
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
                                                 radiusRatioLimit=0.7,
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

    psdEqsp = Measure.getParticleSizeDistribution( psSummary=pss,
                                                    sizeParam='eqsp',
                                                    saveData=True,
                                                    sampleName=dataName,
                                                    outputDir=ofl )

    #-Fabric-
    contactTableRW = Measure.getContactNormalsSPAM(corLabMap, method = 'randomWalker')
    N, F, Fq = Measure.fabricVariablesWithUncertainity( contactTableRW, vectUncert = 0.26 )

    np.savetxt((ofl+ dataName + '-' + str(eLen/d50) +'D50-contactTableRW.csv'), contactTableRW, delimiter=',')    # Contact table RW
    np.savetxt((ofl+ dataName + '-' + 'N.txt'), N, fmt='%r')                                       # Fabric tensor
    np.savetxt((ofl+ dataName + '-' + 'F.txt'), F, fmt='%r')                                       # Deviatoric fabric tensor
    np.savetxt((ofl+ dataName + '-' + 'Fq.txt'), Fq, fmt='%r')                                     # Ansiotropy factor

totalTimeEnd = time.time()
totalTimeTaken = totalTimeEnd - totalTimeStart

print('\n\n--------------------------------------**')
print('Total time taken to analyze(mins): ~' + str(totalTimeTaken//60))
