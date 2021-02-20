'''
Code for the analysis of 3D tomo data for samples collected around the cone penetrometer

Steps:
    -Initialize the samples variables
    -Read the files and extract subregions from the 3D tomo data - 6D50 edge length

    -Compute particle sizes of the sand samples - user threshold after initial guess from OTSU
    -Get particel aspect ration from PCA length ratios
    -Get particle size distribution
    -Calculate relative breakage parameter
    ---
    -Get contact normals
    -Get fabric tensors
    -Get fabric projection plots

'''

from pac import Read
from pac import Filter
from pac import Segment
from pac import Measure
from pac import Plot

v=int(input('Enter number: '))
userH=int(1)
userRR=0.6
mainInput = '/home/eg/codes/pacInput/'
mainOutput = '/home/eg/codes/pacOutput/cone/finalFolder/'

# 2QR 25
if v==1:
    scanList = ['2QR_25_tip','2QR_25_mid','2QR_25_top']
    originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
    psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                    maxRows=10,
                                    dataForm='array').reshape(10,2)
# 2QR 50
if v==2:
    scanList = ['2QR_50_tip','2QR_50_mid','2QR_50_top']
    originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
    psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                    maxRows=10,
                                    dataForm='array').reshape(10,2)

# 2QR 90
if v==3:
    scanList = ['2QR_90_tip','2QR_90_mid','2QR_90_top']
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

numberofSubregionsPerScan=10
nD50=6

for scan in scanList:

    # Locations of data
    scanInputLoc = mainInput + scan + '/'
    subregionInfo = mainInput + scan + '/' + 'subregionInfo.csv'
    outputLoc = mainOutput + scan + '/'

    # Reading subregion data saved in the subregionInfo location.
    subregionCalib = Read.readDataFromCsv( subregionInfo,
                                             skipHeader=1,
                                             maxRows=1,
                                             fmt='float',
                                             dataForm='number' )

    # The average particle size of the sand
    subregionD50 = Read.readDataFromCsv( subregionInfo,
                                           skipHeader=2,
                                           maxRows=1,
                                           fmt='float',
                                           dataForm='number')

    # Slice number for center of the scan
    subregionZ = Read.readDataFromCsv( subregionInfo,
                                         skipHeader=3,
                                         maxRows=1,
                                         fmt='int',
                                         dataForm='number')

    # Location of points in the scan - top left corner of the cube
    subregionYXArray = Read.readDataFromCsv( subregionInfo,
                                               skipHeader=4,
                                               maxRows=numberofSubregionsPerScan,
                                               fmt='int',
                                               dataForm='array' ).reshape(numberofSubregionsPerScan,2)

    for currentSubregion in range(0,numberofSubregionsPerScan):

        currentSampleName = scan + '-' + str( round( currentSubregion ) )

        # Extraction of the subregion from the complete scan
        subregionGLIMap = Read.extractSubregionFromTiffSequence( folderDir=scanInputLoc,
                                                                   reference='topLeft',
                                                                   Z=subregionZ,
                                                                   Y=subregionYXArray[currentSubregion,0],
                                                                   X=subregionYXArray[currentSubregion,1],
                                                                   lngt=nD50*subregionD50,
                                                                   calib=subregionCalib,
                                                                   invImg=False,
                                                                   saveImg=False,
                                                                   outputDir=outputLoc,
                                                                   sampleName=currentSampleName )

        # Binarization using Otsu
        binMap = Segment.binarizeAccordingToOtsu( gliMapToBinarize=subregionGLIMap,
                                                  sampleName=currentSampleName,
                                                  saveImg=False,
                                                  outputDir=outputLoc,
                                                  returnThresholdVal=False)

        # EDM and particle centers
        edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
                                                  scaleUp = int(1),
                                                  saveImg=False,
                                                  sampleName=currentSampleName,
                                                  outputDir=outputLoc )
        # EDM peaks
        edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
                                                        h=userH,
                                                        sampleName=currentSampleName,
                                                        saveImg=False,
                                                        outputDir=outputLoc )

        # Watershed segmentation
        labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
                                                edmMapForTopo=edmMap,
                                                edmPeaksForSeed=edmPeaksMap,
                                                sampleName=currentSampleName,
                                                saveImg=False,
                                                outputDir=outputLoc )

        # Correction of segmentation errors
        corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
                                                     pad=int(2),
                                                     areaLimit=700,
                                                     considerEdgeLabels=True,
                                                     checkForSmallParticles=False,
                                                     radiusCheck=True,
                                                     radiusRatioLimit=userRR,
                                                     sampleName=currentSampleName,
                                                     saveImg=True,
                                                     outputDir=outputLoc )

        # Removal of edge labels
        noEdgeCorLabMap = Segment.removeEdgeLabels( labelledMapForEdgeLabelRemoval=corLabMap,
                                                    pad=0,
                                                    sampleName=currentSampleName,
                                                    saveImg=True,
                                                    outputDir=outputLoc )

        # Particle size list
        pss = Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis=noEdgeCorLabMap,
                                            calibrationFactor=subregionCalib,
                                            saveData=True,
                                            sampleName=currentSampleName,
                                            outputDir=outputLoc )

        # Particle size distribution
        psdFeret = Measure.getParticleSizeDistribution( psSummary=pss,
                                                        sizeParam='feretMin',
                                                        saveData=True,
                                                        sampleName=currentSampleName,
                                                        outputDir=outputLoc )

        # Relative breakage
        brHardin = Measure.getRelativeBreakageHardin( psdOriginal=psdOrig,
                                                      psdCurrent=psdFeret[:,2:],
                                                      smallSizeLimit=0.075,
                                                      saveData=True,
                                                      sampleName=currentSampleName,
                                                      outputDir=outputLoc )
