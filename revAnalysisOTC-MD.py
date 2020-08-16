# Importing files
from pac import Reader                      # Reads files, cuts into smaller pieces
from pac import Filter                      # Filters files using NLM filter
from pac import Segment                     # Binarizatio and WS
from pac import Measure                     # Calculates particle size, mprphology, contact, breakage
from pac import Plot                        # Plotting functions
import time                                 # The fourth dimension
import matplotlib.pyplot as plt             # matplotlib
import skimage.external.tifffile as tf      # scikit-image
import numpy as np                          # numpy

ifl = '/home/eg/codes/pacInput/OTC-MD-0N/'
ofl = '/home/eg/codes/pacOutput/REV/OTC-MD-0N/'
origGSDLoc = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

# Data details 0N:
dataName = 'otc-MD-0N'
measVRSamp = 0.609                                                  # Void ratio measured from 1D compression experiment
d50 = 0.72                                                          # D50 in mm - original gradation
cal = 0.01193                                                       # calibration from CT mm/voxel
zCenter = 513                                                       # Voxel units - center of slice
yCenter = 405                                                       # Voxel units - vertical center
xCenter = 499                                                       # Voxel units - horizontal center
originalGSD = np.loadtxt( origGSDLoc , delimiter=',' )              # Original GSD

analyzeTotalVol = True
analyzeRevSizes = False

if analyzeTotalVol == True :
    eLen = 7.5*d50 # Max edge length in D50s

    # Particle size
    gliMap = Reader.readTiffFileSequence( ifl,zCenter,yCenter,xCenter,eLen,cal,invImg=False )
    binMap, edMap, edPeakMap, labMap = Segment.obtLabMapITKWS( gliMap,measuredVoidRatio=measVRSamp,outputLocation=ofl,edmScaleUp=1,peakEdLimit=5 )
    corLabMap = Segment.fixErrSeg( labMap, pad=2, outputLocation=ofl, areaLimit = 700)
    noEdgeCorLabMap = Segment.removeEdgeLabels( corLabMap )

    # Save files:
    gliNameAll = ofl + 'maxD50-gliMap.tiff'
    binNameAll = ofl + 'maxD50-binMap.tiff'
    edNameAll = ofl + 'maxD50-edMap.tiff'
    labNameAll = ofl + 'maxD50-labMap.tiff'
    corLabNameAll = ofl + 'maxD50-corLabMap.tiff'
    noEdgeCorLabNameAll = ofl + 'maxD50-noEdgeCorLabMap.tiff'
    tf.imsave( gliNameAll,gliMap.astype( 'uint32' ) )
    tf.imsave( binNameAll,binMap.astype( 'uint32' ) )
    tf.imsave( edNameAll, edMap.astype( 'uint32' ) )
    tf.imsave( labNameAll,labMap.astype( 'uint32' ) )
    tf.imsave( corLabNameAll, corLabMap.astype( 'uint32' ) )
    tf.imsave( noEdgeCorLabNameAll, noEdgeCorLabMap.astype( 'uint32' ) )

    # Grain size distribution
    gsd1a, gsd2a, gas3a, gas4a, gas5a, gas6a = Measure.gsdAll( noEdgeCorLabMap, calib = cal )
    np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd1.csv'), gsd1a, delimiter=',')
    np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd2.csv'), gsd2a, delimiter=',')
    np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd3.csv'), gsd3a, delimiter=',')
    np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd4.csv'), gsd4a, delimiter=',')
    np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd5.csv'), gsd5a, delimiter=',')
    np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd6.csv'), gsd6a, delimiter=',')
    Plot.grainSizeDistribution( originalGSD, gsd1a, gsd2a, gas3a, gas4a, gas5a, gas6a )

    # Contact and fabric
    contactTableRWAll = Measure.contactNormalsSpam( corLabMap, method='rw')
    np.savetxt( ( ofl + 'maxD50-contactTableRW.csv'), contactTableRWAll, delimiter=',')    # Contact table
    N, F, Fq = Measure.fabricVariablesWithUncertainity( contactTableRWAll, vectUncert = 0.26 )
    np.savetxt( (ofl + 'maxD50-N.txt' ), N, fmt = '%r')    # Fabric tensor
    np.savetxt( (ofl + 'maxD50-F.txt' ), F, fmt = '%r')    # Deviatoric fabric tensor
    np.savetxt( (ofl + 'maxD50-Fq.txt' ), Fq, fmt = '%r')  # Ansiotropy factor

if analyzeRevSizes == True :
    if analyzeTotalVol == False:
        binMap = tf.imread(binNameAll).astype('uint32')
        gliMap = tf.imread(gliNameAll).astype('uint32')

    sizeRange = np.arange(2, 8, 1)
    zCenterSu = binMap.shape[0]//2
    yCenterSu = binMap.shape[1]//2
    xCenterSu = binMap.shape[2]//2

    sizeList = []
    voidRatioList = []
    gsdList = []

    for i in sizeRange:
        length = i*d50
        sizeList.append(i)

        print('\nAnalyzing subregion of size : ' + str(i) + 'd50')

        upSlice = int(zCenterSu + round( ( length / 2 ) / cal ))
        lowSlice = int(zCenterSu - round( ( length / 2 ) / cal ))
        upRow = int(yCenterSu + round( ( length / 2 ) / cal ))
        lowRow = int(yCenterSu - round( ( length / 2 ) / cal ))
        upCol = int(xCenterSu + round( ( length / 2 ) / cal ))
        lowCol = int(xCenterSu - round( ( length / 2 ) / cal ))

        subRegionBinMap = binMap[lowSlice:upSlice,lowRow:upRow,lowCol:upCol]
        subVoidRatio = Segment.calcVoidRatio(subRegionBinMap)
        voidRatioList.append(subVoidRatio)

        subRegionGliMap = gliMap[lowSlice:upSlice,lowRow:upRow,lowCol:upCol]
        binMask, edMap, edPeaksMap, subRegionLabMap = Segment.obtLabMapITKWS( subRegionGliMap,
                                                                              knownThreshold=10574.2,
                                                                              outputLocation=ofl,
                                                                              edmScaleUp=1,
                                                                              peakEdLimit=5)

        subRegionCorrectedLabMap = Segment.fixErrSeg( subRegionLabMap,
                                                      pad=2,
                                                      outputLocation=ofl,
                                                      areaLimit=700 )

        subRegionNoEdgeCorrectedLabMap = Segment.removeEdgeLabels( subRegionCorrectedLabMap )

        gsd1,gsd2,gsd3,gsd4,gsd5,gsd6 = Measure.gsdAll( subRegionNoEdgeCorrectedLabMap , calib=cal )
        np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd1.csv'), gsd1, delimiter=',')
        np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd2.csv'), gsd2, delimiter=',')
        np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd3.csv'), gsd3, delimiter=',')
        np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd4.csv'), gsd4, delimiter=',')
        np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd5.csv'), gsd5, delimiter=',')
        np.savetxt( ( ofl + str( np.round( i ) ) +'D50-gsd6.csv'), gsd6, delimiter=',')

        contactTableRW = Measure.contactNormalsSpam(subRegionCorrectedLabMap, method='rw')
        np.savetxt( ( ofl + str( np.round( i ) ) +'D50-contactTableRW.csv'), contactTableRW, delimiter=',')    # Contact table RW
        N, F, Fq = Measure.fabricVariablesWithUncertainity( contactTableRW, vectUncert = 0.26 )
        np.savetxt((ofl + str( np.round( i ) ) +'-D50-N.txt'), N, fmt='%r')    # Fabric tensor
        np.savetxt((ofl + str( np.round( i ) ) +'-D50-F.txt'), F, fmt='%r')    # Deviatoric fabric tensor
        np.savetxt((ofl + str( np.round( i ) ) +'-D50-Fq.txt'), Fq, fmt='%r')  # Ansiotropy factor

    # Plot grain size distribution
    plt.plot(sizeList,voidRatioList)
    plt.ylim([0,2])
    plt.xlim([0,10])
    plt.show()

    # Save void ratio list as output
    textFileName = ofl + 'voidRatio.txt'
    with open(textFileName, "w") as output:
        output.write(str(voidRatioList))

