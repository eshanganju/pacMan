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

scan = '2QR_25_tip'

numberofSubregionsPerScan=10
nD50=6

mainInput = '/home/eg/codes/pacInput/'
mainOutput = '/home/eg/codes/pacOutput/cone/'

# Locations of data
scanInputLoc = mainInput + scan + '/'
subregionInfo = mainInput + scan + '/' + 'subregionInfo.csv'
outputLoc = mainOutput + scan + '/'

# Reading subregion data saved in the subregionInfo location.
subregionCalib = Reader.readDataFromCsv( subregionInfo,
                                         skipHeader=1,
                                         skipFooter=numberofSubregionsPerScan+2,
                                         fmt='float',
                                         dataForm='number' )

# The average particle size of the sand
subregionD50 = Reader.readDataFromCsv( subregionInfo,
                                       skipHeader=2,
                                       skipFooter=numberofSubregionsPerScan+1,
                                       fmt='float',
                                       dataForm='number')

# Slice number for center of the scan
subregionZ = Reader.readDataFromCsv( subregionInfo,
                                     skipHeader=3,
                                     skipFooter=numberofSubregionsPerScan,
                                     fmt='int',
                                     dataForm='number')

# Location of points in the scan - top left corner of the cube
subregionYXArray = Reader.readDataFromCsv( subregionInfo,
                                           skipHeader=4,
                                           fmt='int',
                                           dataForm='array' )

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
                                                    h=int(5),                       # 5 works well
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
                                                 radiusRatioLimit=0.7,
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
    psList = Segment.getParticleSize( labelledMapForParticleSizeAnalysis=noEdgeCorLabMap,
                                      calibrationFactor=subregionCalib,
                                      sampleName=currentSampleName,
                                      saveData=True,
                                      outputDir=outputLoc)


