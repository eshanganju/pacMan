"""Analysis of  large amount of crushing
"""

from pac import Reader
from pac import Segment
from pac import Measure
from pac import Plot

v = int( input( 'enter number' ) )

if v==0:
    centerZ=513
    centerY=488
    centerX=493
    d50=0.73
    calibVal=0.01193
    sampleName='2QR_D_90'
    inputDir='/home/eg/codes/pacInput/2QR-4500N/'
    outputDir='/home/eg/codes/pacOutput/1D/90MPaAnalysis/2QR_D_90/'
    originalPSDLoc='/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
    originalPSD=Reader.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

if v==1:
    centerZ=508
    centerY=505
    centerX=489
    d50=0.62
    calibVal=0.01193
    sampleName='OGF_D_90'
    inputDir='/home/eg/codes/pacInput/OGF-4500N/'
    outputDir='/home/eg/codes/pacOutput/1D/90MPaAnalysis/OGF_D_90/'
    originalPSDLoc='/home/eg/codes/pacInput/originalGSD/ogfOrig.csv'
    originalPSD=Reader.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

if v==2:
    centerZ=513
    centerY=497
    centerX=502
    d50=0.72
    calibVal=0.01193
    sampleName='OTC_D_90'
    inputDir='/home/eg/codes/pacInput/OTC-4500N/'
    outputDir='/home/eg/codes/pacOutput/1D/90MPaAnalysis/OTC_D_90/'
    originalPSDLoc='/home/eg/codes/pacInput/originalGSD/otcOrig.csv'
    originalPSD=Reader.readDataFromCsv( originalPSDLoc,
                                        maxRows=10,
                                        dataForm='array').reshape(10,2)

gliMap = Reader.extractSubregionFromTiffSequence( folderDir=inputDir,
                                                  reference='center',
                                                  Z=centerZ,
                                                  Y=centerY,
                                                  X=centerX,
                                                  lngt=6*d50,
                                                  calib=calibVal,
                                                  invImg=False,
                                                  outputDir=outputDir,
                                                  sampleName=sampleName )

binMap = Segment.binarizeAccordingToOtsu( gliMapToBinarize=gliMap,
                                          sampleName=sampleName,
                                          saveImg=False,
                                          outputDir=outputDir,
                                          returnThresholdVal=False )

edtMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=binMap,
                                          scaleUp=int(1),
                                          saveImg=False,
                                          sampleName=sampleName,
                                          outputDir=outputDir )

edpMap = Segment.obtainLocalMaximaMarkers( edtMap,
                                           h=int(1),
                                           sampleName=sampleName,
                                           saveImg=False,
                                           outputDir=outputDir )

labMap = Segment.segmentUsingWatershed( binaryMapToSeg=binMap,
                                        edmMapForTopo=edtMap,
                                        edmPeaksForSeed=edpMap,
                                        sampleName=sampleName,
                                        saveImg=False,
                                        outputDir=outputDir )

corLabMap = Segment.fixErrorsInSegmentation( labMap,
                                             pad=int(2),
                                             areaLimit=700,
                                             considerEdgeLabels=True,
                                             checkForSmallParticles=False,
                                             radiusCheck=True,
                                             radiusRatioLimit=0.6,
                                             sampleName=sampleName,
                                             saveImg=True,
                                             outputDir=outputDir )

noEdgeCorLabMap = Segment.removeEdgeLabels( corLabMap,
                                            pad=0,
                                            sampleName=sampleName,
                                            saveImg=True,
                                            outputDir=outputDir )

pss = Measure.getParticleSizeArray( noEdgeCorLabMap,
                                    calibrationFactor=calibVal,
                                    saveData=True,
                                    sampleName=sampleName,
                                    outputDir=outputDir )

psdFeret = Measure.getParticleSizeDistribution( psSummary=pss,
                                                sizeParam='feretMin',
                                                saveData=True,
                                                sampleName=sampleName,
                                                outputDir=outputDir )

brHardin = Measure.getRelativeBreakageHardin( psdOriginal=originalPSD,
                                              psdCurrent=psdFeret[:,2:],
                                              smallSizeLimit=0.075,
                                              saveData=True,
                                              sampleName=sampleName,
                                              outputDir=outputDir )

