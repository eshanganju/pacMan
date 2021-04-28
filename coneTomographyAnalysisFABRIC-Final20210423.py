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

import numpy as np

import tifffile as tf

userH=int(1)
userRR=0.8
nD50=6

mainInput = '/home/eg/codes/pacInput/resinSampleUpdate/'
subregionInput = '/home/eg/pacMan/subregionInfoResin/'
mainOutput = '/home/eg/codes/pacOutput/cone/finalFolderResin/'

originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv'
psdOrig = Read.readDataFromCsv( originalPSDLoc,maxRows=10,dataForm='array').reshape(10,2)

# scanList = [1,5]
scanList = [2,3,4]

for scan in scanList:

    scanInputLoc = mainInput + 'scan' + str(int(scan)) + '/corData/'
    subregionInfo = subregionInput + 'subregionInfo-' + 'scan' + str(int(scan)) +'.csv'
    outputLoc = mainOutput + 'scan' + str(int(scan))  + '/'

    if scan == 1: numberofSubregionsPerScan=9
    if scan == 2: numberofSubregionsPerScan=6
    if scan == 3: numberofSubregionsPerScan=6
    if scan == 4: numberofSubregionsPerScan=6
    if scan == 5: numberofSubregionsPerScan=9

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
                                         dataForm='number' )

    # Slice number for center of the scan
    subregionZ = Read.readDataFromCsv( subregionInfo,
                                         skipHeader=3,
                                         maxRows=1,
                                         fmt='int',
                                         dataForm='number' )

    # Location of points in the scan - top left corner of the cube
    subregionYXArray = Read.readDataFromCsv( subregionInfo,
                                               skipHeader=4,
                                               maxRows=numberofSubregionsPerScan,
                                               fmt='int',
                                               dataForm='array' ).reshape( numberofSubregionsPerScan,2 ) 

    start=0
    for currentSubregion in range( start, numberofSubregionsPerScan ):

        currentSampleName = 'scan' + str(int(scan)) + '-' + str( round( currentSubregion ) )

        # Extraction of the subregion from the complete scan
        subregionGLIMap = Read.extractSubregionFromTiffSequence( folderDir=scanInputLoc,
                                                                 reference='topLeft',
                                                                 Z=subregionZ,
                                                                 Y=subregionYXArray[currentSubregion,0],
                                                                 X=subregionYXArray[currentSubregion,1],
                                                                 lngt=nD50*subregionD50,
                                                                 calib=subregionCalib,
                                                                 invImg=False,
                                                                 saveImg=True,
                                                                 outputDir=outputLoc,
                                                                 sampleName=currentSampleName )
        
        # Filter using NLM
        filteredSubregionGLIMap = Filter.filterUsingNlm( gli=subregionGLIMap,
                                                         bitDepth=16,
                                                         pSize=5,
                                                         pDistance=7,
                                                         saveImg=True,
                                                         outputDir=outputLoc,
                                                         sampleName=currentSampleName )

        # Binarization using Otsu
        binMap = Segment.binarizeAccordingToOtsu( gliMapToBinarize=filteredSubregionGLIMap,
                                                  sampleName=currentSampleName,
                                                  saveImg=True,
                                                  outputDir=outputLoc,
                                                  returnThresholdVal=False )

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
                                                     checkForSmallParticles=True, # This has also been changed since the last code
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

        # Relative breakage - this section will have to be updated
        brHardin = Measure.getRelativeBreakageHardin( psdOriginal=psdOrig,
                                                      psdCurrent=psdFeret[:,2:],
                                                      smallSizeLimit=0.075,
                                                      saveData=True,
                                                      sampleName=currentSampleName,
                                                      outputDir=outputLoc )

        # Contact analysis - With edge labels included
        '''The "corlabMap" is read again to prevent errors in stored memory.
        '''
        corLabMapLoc = outputLoc + currentSampleName + '-correctedLabelMap.tif'
        corLabMap = tf.imread(corLabMapLoc)
            
        coordinationNumberList = Measure.getCoordinationNumberList( corLabMap, excludeEdgeLabels=True )
        np.savetxt( ( outputLoc + currentSampleName + '-CNList(edgeLabelsExcluded).csv'), coordinationNumberList, delimiter=',')

        contactTableRW = Measure.getContactNormalsSPAM(corLabMap, method = 'randomWalker')
        np.savetxt( ( outputLoc + currentSampleName +'-contactTable-RW.csv'), contactTableRW, fmt='%r' )     # Contact table
        
        N, F, Fq = Measure.fabricVariablesWithUncertainity( contactTableRW, vectUncert = 0.26 )
        np.savetxt( ( outputLoc + currentSampleName +'-N2.txt'), N, fmt='%r' )                                   # Fabric tensor
        np.savetxt( ( outputLoc + currentSampleName +'-F2.txt'), F, fmt='%r' )                                   # Deviatoric fabric tensor
        np.savetxt( ( outputLoc + currentSampleName +'-Fq2.txt'), Fq, fmt='%r' )                                 # Ansiotropy factor



