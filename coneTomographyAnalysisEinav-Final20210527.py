'''
This code is to calculate the particle size distribution along with removal of particles smaller than threshold for previous analysis
theFollowing thisg will be done: 
    1. Read the noEdgeCorLabMap
    2. remove small particles from the noEdgeCorLabMap (measure loss + )
    3. Carry out particle size analysis

'''

from pac import Read
from pac import Filter
from pac import Segment
from pac import Measure
from pac import Plot

import tifffile as tf

userH=int(1)
userRR=0.8 # This has been updated from the number used previously (0.6). Check the psds
# numberofSubregionsPerScan=10
nD50=6

mainInput = '/home/eg/codes/pacInput/'
subregionInput = '/home/eg/pacMan/subregionInfo/'
mainOutput = '/home/eg/codes/pacOutput/cone/finalFolder/'

for v in [5]:
    # 2QR 25
    if v==1:
        scanList = ['2QR_25_tip','2QR_25_mid','2QR_25_top']
        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)
    # 2QR 50
    if v==2:
        scanList = ['2QR_50_tip']

        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

    # 2QR 90
    if v==3:
        scanList = ['2QR_90_tip']

        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

    # OGF 25
    if v==4:
        scanList = ['OGF_25_tip','OGF_25_mid','OGF_25_top']
        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

    # OGF 50
    if v==5:
        scanList = ['OGF_50_tip','OGF_50_mid','OGF_50_top']
        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

    # OGF 90
    if v==6:
        scanList = ['OGF_90_tip','OGF_90_mid','OGF_90_top']
        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

    # OTC 25
    if v==7:
        scanList = ['OTC_25_tip','OTC_25_mid','OTC_25_top']
        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)
    # OTC 50
    if v==8:
        scanList = ['OTC_50_tip','OTC_50_mid','OTC_50_top']
        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

    # OTC 90
    if v==9:
        scanList = ['OTC_90_tip','OTC_90_mid','OTC_90_top']
        originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv'
        psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

    for scan in scanList:

        # Locations of data
        scanInputLoc = mainInput + scan + '/'
        subregionInfo = subregionInput + 'subregionInfo-' + scan +'.csv'
        outputLoc = mainOutput + scan + '/'

        # Reading subregion data saved in the subregionInfo location.
        subregionCalib = Read.readDataFromCsv( subregionInfo,
                                                 skipHeader=1,
                                                 maxRows=1,
                                                 fmt='float',
                                                 dataForm='number' )

        # # The average particle size of the sand
        # subregionD50 = Read.readDataFromCsv( subregionInfo,
        #                                        skipHeader=2,
        #                                        maxRows=1,
        #                                        fmt='float',
        #                                        dataForm='number')

        # # Slice number for center of the scan
        # subregionZ = Read.readDataFromCsv( subregionInfo,
        #                                      skipHeader=3,
        #                                      maxRows=1,
        #                                      fmt='int',
        #                                      dataForm='number')

        # # Location of points in the scan - top left corner of the cube
        # subregionYXArray = Read.readDataFromCsv( subregionInfo,
        #                                            skipHeader=4,
        #                                            maxRows=numberofSubregionsPerScan,
        #                                            fmt='int',
        #                                            dataForm='array' ).reshape(numberofSubregionsPerScan,2)

        # 2QR 90 tip was not run completely.
        if scan == '2QR_90_tip': numberofSubregionsPerScan = 7
        else: numberofSubregionsPerScan = 10
        
        for currentSubregion in range(0,numberofSubregionsPerScan):

            currentSampleName = scan + '-' + str( round( currentSubregion ) )

            # # Extraction of the subregion from the complete scan
            # subregionGLIMap = Read.extractSubregionFromTiffSequence( folderDir=scanInputLoc,
            #                                                            reference='topLeft',
            #                                                            Z=subregionZ,
            #                                                            Y=subregionYXArray[currentSubregion,0],
            #                                                            X=subregionYXArray[currentSubregion,1],
            #                                                            lngt=nD50*subregionD50,
            #                                                            calib=subregionCalib,
            #                                                            invImg=False,
            #                                                            saveImg=False,
            #                                                            outputDir=outputLoc,
            #                                                            sampleName=currentSampleName )

            # # Binarization using Otsu
            # binMap = Segment.binarizeAccordingToOtsu( gliMapToBinarize=subregionGLIMap,
            #                                           sampleName=currentSampleName,
            #                                           saveImg=False,
            #                                           outputDir=outputLoc,
            #                                           returnThresholdVal=False)

            # # EDM and particle centers
            # edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
            #                                           scaleUp = int(1),
            #                                           saveImg=False,
            #                                           sampleName=currentSampleName,
            #                                           outputDir=outputLoc )
            # # EDM peaks
            # edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
            #                                                 h=userH,
            #                                                 sampleName=currentSampleName,
            #                                                 saveImg=False,
            #                                                 outputDir=outputLoc )

            # # Watershed segmentation
            # labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
            #                                         edmMapForTopo=edmMap,
            #                                         edmPeaksForSeed=edmPeaksMap,
            #                                         sampleName=currentSampleName,
            #                                         saveImg=False,
            #                                         outputDir=outputLoc )

            # # Correction of segmentation errors
            # corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
            #                                              pad=int(2),
            #                                              areaLimit=700,
            #                                              considerEdgeLabels=True,
            #                                              checkForSmallParticles=True, # This has also been changed since the last code
            #                                              radiusCheck=True,
            #                                              radiusRatioLimit=userRR,
            #                                              sampleName=currentSampleName,
            #                                              saveImg=True,
            #                                              outputDir=outputLoc )

            # # Removal of edge labels
            # noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
            #                                             pad=0,
            #                                             sampleName=currentSampleName,
            #                                             saveImg=True,
            #                                             outputDir=outputLoc )

            noEdgeCorlabMap = tf.imread( outputLoc + currentSampleName + '-noEdgeCorrectedLabelMap.tif' )

            smallPtclRemNoEdgeCorlabMap = Segment.removeSmallParticles( labMapWithSmallPtcl=noEdgeCorlabMap,
                                                                        voxelCountThreshold = 1000, 
                                                                        saveImg=True, 
                                                                        sampleName=currentSampleName, 
                                                                        outputDir=outputLoc)

            noEdgeCorlabMap = tf.imread( outputLoc + currentSampleName + '-noEdgeCorrectedLabelMap.tif' )

            solidVolumeWithSmallParticles, solidVolumeWithoutSmallParticles, volumeLoss, percentLoss = Segment.computeVolumeLoss( labelledMapWithSmallParticles = noEdgeCorlabMap, 
                                                                                                                                  labelledMapWithoutSmallParticles = smallPtclRemNoEdgeCorlabMap, 
                                                                                                                                  saveFile=True, 
                                                                                                                                  sampleName=currentSampleName, 
                                                                                                                                  outputDir=outputLoc )

            # Particle size list
            pss = Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis=smallPtclRemNoEdgeCorlabMap,
                                                calibrationFactor=subregionCalib,
                                                saveData=True,
                                                sampleName=currentSampleName + '-noSmallPtclSPAMFeret',
                                                outputDir=outputLoc )

            # Particle size distribution
            psdFeret = Measure.getParticleSizeDistribution( psSummary=pss,
                                                            sizeParam='feretMin',
                                                            saveData=True,
                                                            sampleName=currentSampleName + '-noSmallPtclSPAMFeret',
                                                            outputDir=outputLoc )

            # # Relative breakage - this section will have to be updated
            # brHardin = Measure.getRelativeBreakageHardin( psdOriginal=psdOrig,
            #                                               psdCurrent=psdFeret[:,2:],
            #                                               smallSizeLimit=0.075,
            #                                               saveData=True,
            #                                               sampleName=currentSampleName + '-noSmallPtclSPAMFeret',
            #                                               outputDir=outputLoc )
