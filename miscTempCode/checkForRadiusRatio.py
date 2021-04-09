"""Checking values of radius ratio
        The h value is set at 1 in the segmentation
        the rr is set by user.

"""

import tifffile as tiffy

from pac import Segment
from pac import Measure
from pac import Read

userRR=float(input('Enter rr: '))
currentSampleName='test-1-' + str(userRR)
outputLoc='/home/eg/codes/pacOutput/cone/testForOversegmentationCorrection/testOfPeaksAndRr/'

labMap = tiffy.imread('/home/eg/codes/pacOutput/cone/testForOversegmentationCorrection/testOfPeaksAndRr/test-1-0.6-labMap.tif')

subregionCalib = 0.0123007
originalPSDLoc = '/home/eg/codes/pacInput/originalGSD/2qrOrig.csv'
psdOrig = Read.readDataFromCsv( originalPSDLoc,
                                maxRows=10,
                                dataForm='array').reshape(10,2)

corLabMap = Segment.fixErrorsInSegmentation( labelledMapForOSCorr=labMap,
                                             pad=int(2),
                                             areaLimit=700,
                                             considerEdgeLabels=True,
                                             checkForSmallParticles=True,
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

