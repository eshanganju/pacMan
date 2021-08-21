'''
This code is for the analysis offabric in an air pluviated OTC sand sample

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

mainInput = '/home/eg/codes/pacInput/AP-OTC-85_DR-HC/'
subregionInput = '/home/eg/pacMan/subregionInfoAP-OTC/'
mainOutput = '/home/eg/codes/pacOutput/HC/'

originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv'
psdOrig = Read.readDataFromCsv( originalPSDLoc,maxRows=10,dataForm='array').reshape(10,2)

scanList = [1]

for scan in scanList:

    scanInputLoc = mainInput + 'images/'
    subregionInfo = subregionInput + 'subregionInfo-AP-OTC.csv'
    outputLoc = mainOutput + 'AP-D-OTC/'

    if scan == 1: numberofSubregionsPerScan=1

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

    if True:

       	currentSubregion = 0
       	currentSampleName = 'AP-D-OTC'

        print('\n\nRunning: ' + currentSampleName)

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
        
        # # Filter using NLM
        # filteredSubregionGLIMap = Filter.filterUsingNlm( gli=subregionGLIMap,
        #                                                  bitDepth=16,
        #                                                  pSize=5,
        #                                                  pDistance=7,
        #                                                  saveImg=True,
        #                                                  outputDir=outputLoc,
        #                                                  sampleName=currentSampleName,
        #                                                  loop=False )

        # Binarization using Otsu
        binMap = Segment.binarizeAccordingToOtsu( gliMapToBinarize=subregionGLIMap,
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
                                                saveImg=True,
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