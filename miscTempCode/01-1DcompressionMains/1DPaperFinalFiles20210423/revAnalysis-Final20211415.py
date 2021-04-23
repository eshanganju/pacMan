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

timeStart = time.time()

run = ['OGF-0N', 'OTC-MD-0N','2QR-0N-Middle', 'OTC-0N' ]

for i in run:

    # Void ratios taken from /1DAnalysisFinal20210415/VoidRatioCalculationsCTSamples.ods

    ############# 2QR-D #############   
    if i == '2QR-0N-Middle' :
        inputFolderLocation = '/home/eg/codes/pacInput/2QR-0N/'
        ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/REV/2QR-0N-Middle/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv' # Original GSD location

        dataName = i
        measVoidRatio = 0.734                                     # Void ratio measured from 1D compression experiment
        d50 = 0.73                                                          # D50 in mm - original gradation
        cal = 0.011932                                                      # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 427                                                       # Voxel units - vertical center
        xCenter = 504                                                       # Voxel units - horizontal center
        origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )         # Original GSD
        mainMapnD50 = 7.5
        eLen = mainMapnD50*d50

    if i == 'OGF-0N' :
        inputFolderLocation = '/home/eg/codes/pacInput/OGF-0N/'
        ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/REV/OGF-0N/' # output folder location ofl
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
        mainMapnD50 = 8.5
        eLen = mainMapnD50*d50

    if i == 'OTC-0N' :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-0N/'
        ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/REV/OTC-0N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details 0N:
        dataName = i
        measVoidRatio = 0.534                                               # Void ratio measured from 1D compression experiment
        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 457                                                       # Voxel units - vertical center
        xCenter = 507                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )          # Original GSD
        mainMapnD50 = 7.5
        eLen = mainMapnD50*d50

    if i == 'OTC-MD-0N' :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-0N/'
        ofl = '/home/eg/codes/pacOutput/1DAnalysisFinal20210415/REV/OTC-MD-0N/'
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
        mainMapnD50 = 7.5
        eLen = mainMapnD50*d50

    mainGliMap = Read.extractSubregionFromTiffSequence( inputFolderLocation,
                                                        Z=zCenter,
                                                        Y=yCenter,
                                                        X=xCenter,
                                                        lngt=eLen,
                                                        calib=cal,
                                                        reference='center',
                                                        invImg=False,
                                                        saveImg=True, 
                                                        outputDir=ofl, 
                                                        sampleName=(dataName+'-Main'))

    binThreshold, mainBinMap = Segment.binarizeAccordingToDensity( gliMapToBinarize=mainGliMap,
                                                                   measuredVoidRatio=measVoidRatio,
                                                                   returnThresholdVal=True,
                                                                   saveImg=True,
                                                                   saveData=True,
                                                                   sampleName=(dataName+'-Main'),
                                                                   outputDir=ofl)

    for nD50 in [7, 6, 5, 4, 3, 2 , 1]:
        eLen = d50 * nD50
        subDataName = dataName + '-' + str(round(nD50)) + 'nD50'
        
        gliMap = Read.extractSubregionFromTiffSequence( inputFolderLocation,
                                                        Z=zCenter,
                                                        Y=yCenter,
                                                        X=xCenter,
                                                        lngt=eLen,
                                                        calib=cal,
                                                        reference='center',
                                                        invImg=False,
                                                        saveImg=True, 
                                                        outputDir=ofl, 
                                                        sampleName=subDataName) 

        binMap = Segment.binarizeAccordingToUserThreshold( gliMapToBinarize=gliMap,
                                                           userThreshold=binThreshold, 
                                                           returnThresholdVal=False, 
                                                           saveImg=True, 
                                                           saveData=True, 
                                                           sampleName=subDataName, 
                                                           outputDir=ofl)

        voidRatio = Segment.calcVoidRatio(binMap)
        voidRatioFile = open( ofl + subDataName + '-voidRatio.Txt',"a+")
        voidRatioFile.write(str(voidRatio))
        voidRatioFile.close()

        edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
                                                  scaleUp = int(1),
                                                  saveImg=False,
                                                  sampleName=subDataName,
                                                  outputDir=ofl )

        # EDM peaks
        edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
                                                        h=1,
                                                        sampleName=subDataName,
                                                        saveImg=False,
                                                        outputDir=ofl )

        # Watershed segmentation
        labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
                                                edmMapForTopo=edmMap,
                                                edmPeaksForSeed=edmPeaksMap,
                                                sampleName=subDataName,
                                                saveImg=True,
                                                outputDir=ofl )

        # Correction of segmentation errors
        corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
                                                     pad=int(2),
                                                     areaLimit=700,
                                                     considerEdgeLabels=True,
                                                     checkForSmallParticles=True,
                                                     radiusCheck=True,
                                                     radiusRatioLimit=0.8,
                                                     sampleName=subDataName,
                                                     saveImg=True,
                                                     outputDir=ofl )

        # Removal of edge labels
        noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
                                                    pad=0,
                                                    sampleName=subDataName,
                                                    saveImg=True,
                                                    outputDir=ofl )

        # Particle size list
        pss = Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis=noEdgeCorLabMap,
                                            calibrationFactor=cal,
                                            saveData=True,
                                            sampleName=subDataName,
                                            outputDir=ofl )

        # Particle size distribution
        psdFeret = Measure.getParticleSizeDistribution( psSummary=pss,
                                                        sizeParam='feretMin',
                                                        saveData=True,
                                                        sampleName=subDataName,
                                                        outputDir=ofl )

timeEnd = time.time()
print('Time taken: ' + str((timeEnd-timeStart)//60))
