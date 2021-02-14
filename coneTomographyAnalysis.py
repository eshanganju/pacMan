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

from pac import Reader
from pac import Filter
from pac import Segment
from pac import Measure
from pac import Plot

v=0

if v==0:
    scanList = ['2QR_25_mid']
    originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
    psdOrig = Reader.readDataFromCsv( originalPSDLoc,
                                      maxRows=10,
                                      dataForm='array').reshape(10,2)

numberofSubregionsPerScan=10
nD50=6

for scan in scanList:
    mainInput = '/home/eg/codes/pacInput/'
    mainOutput = '/home/eg/codes/pacOutput/cone/testEntireNumba/'


    # Locations of data
    scanInputLoc = mainInput + scan + '/'
    subregionInfo = mainInput + scan + '/' + 'subregionInfo.csv'
    outputLoc = mainOutput + scan + '/'

    # Reading subregion data saved in the subregionInfo location.
    subregionCalib = Reader.readDataFromCsv( subregionInfo,
                                             skipHeader=1,
                                             maxRows=1,
                                             fmt='float',
                                             dataForm='number' )

    # The average particle size of the sand
    subregionD50 = Reader.readDataFromCsv( subregionInfo,
                                           skipHeader=2,
                                           maxRows=1,
                                           fmt='float',
                                           dataForm='number')

    # Slice number for center of the scan
    subregionZ = Reader.readDataFromCsv( subregionInfo,
                                         skipHeader=3,
                                         maxRows=1,
                                         fmt='int',
                                         dataForm='number')

    # Location of points in the scan - top left corner of the cube
    subregionYXArray = Reader.readDataFromCsv( subregionInfo,
                                               skipHeader=4,
                                               maxRows=numberofSubregionsPerScan,
                                               fmt='int',
                                               dataForm='array' ).reshape(numberofSubregionsPerScan,2)

    for currentSubregion in range(0,numberofSubregionsPerScan):

        currentSampleName = scan + '-' + str( round( currentSubregion ) )

        # Extraction of the subregion from the complete scan
        subregionGLIMap = Reader.extractSubregionFromTiffSequence( folderDir=scanInputLoc,
                                                                   centerZ=subregionZ,
                                                                   topLeftY=subregionYXArray[currentSubregion,0],
                                                                   topLeftX=subregionYXArray[currentSubregion,1],
                                                                   lngt=nD50*subregionD50,
                                                                   calib=subregionCalib,
                                                                   invImg=False,
                                                                   saveImg=True,
                                                                   outputDir=outputLoc,
                                                                   sampleName=currentSampleName )

        # Binarization using Otsu
        binMap = Segment.binarizeAccordingToOtsu( gliMapToBinarize=subregionGLIMap,
                                                  sampleName=currentSampleName,
                                                  saveImg=True,
                                                  outputDir=outputLoc,
                                                  returnThresholdVal=False)

        # EDM and particle centers
        edmMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
                                                  scaleUp = int(1),
                                                  saveImg=True,
                                                  sampleName=currentSampleName,
                                                  outputDir=outputLoc )
        # EDM peaks
        edmPeaksMap = Segment.obtainLocalMaximaMarkers( edMapForPeaks=edmMap,
                                                        h=int(1),
                                                        sampleName=currentSampleName,
                                                        saveImg=True,
                                                        outputDir=outputLoc )

        # Watershed segmentation
        labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
                                                edmMapForTopo=edmMap,
                                                edmPeaksForSeed=edmPeaksMap,
                                                sampleName=currentSampleName,
                                                saveImg=True,
                                                outputDir=outputLoc )

        # Correction of segmentation errors
        corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
                                                     pad=int(2),
                                                     areaLimit=700,
                                                     considerEdgeLabels=True,
                                                     checkForSmallParticles=True,
                                                     radiusCheck=True,
                                                     radiusRatioLimit=0.6,
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
                                                      psdCurrent=psdFeret,
                                                      smallSizeLimit=0.075,
                                                      saveData=True,
                                                      sampleName=currentSampleName,
                                                      outputDir=outpurLoc )
