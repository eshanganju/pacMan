'''
The objective of this code is to assemble the codes in pac directory to:
    1. Read the data files for the deben scans and extract volume
    2. Segment the data (binarize and WS)
    3. Analyze CT data for
        3.1 Correct REV size
        3.2 Particle size for REV size
        3.3 Particle morphology
        3.4 Relative breakage
        3.5 Fabric tensor
    4. Distribution of rel. breakage and fabric tensor
'''
# Importing files
import sys
sys.path.insert(1, '/home/eg/codes/pac')
from pac import Reader                      	# Reads files, cuts into smaller pieces
from pac import Filter                      	# Filters files using NLM filter
from pac import Segment                     	# Binarizatio and WS
from pac import Measure                     	# Calculates particle size, mprphology, contact, breakage
from pac import Plot                        	# Plotting functions

import time                                						# The fourth dimension
import matplotlib.pyplot as plt            					 	# matplotlib
import skimage.external.tifffile as tf      					# scikit-image
import numpy as np                          					# numpy
from uncertainties import unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *

'''
Binarization according to density measurement
Segmentation according to ITK and auto correction
Particle size values (4) and appropriate gradation
Relative breakage according to Einav
Contact according to ITK and RW
Plotting orientations in rose and EAP diagrams
'''

totalTimeStart=time.time()

# 0 (0 MPa), 500 (10 MPa), 1500 (30 MPa), 4500 (90 MPa)
data = np.array([4500])

for i in data:
    if i == 0 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-0N/'
        ofl = '/home/eg/codes/pacOutput/OTC-0N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details 0N:
        dataName = 'otc-0N'
        measVoidRatio = 0.541                                           # Void ratio measured from 1D compression experiment
        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 457                                                       # Voxel units - vertical center
        xCenter = 507                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

    elif i == 500 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-500N/'
        ofl = '/home/eg/codes/pacOutput/OTC-500N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details:
        dataName = 'otc-500N'
        measVoidRatio = 0.517                                           # Void ratio measured from 1D compression experiment

        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 455                                                       # Voxel units - vertical center
        xCenter = 507                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

    elif i == 1500 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-1500N/'
        ofl = '/home/eg/codes/pacOutput/OTC-1500N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details:
        dataName = 'otc-1500N'
        measVoidRatio = 0.499                                           # Void ratio measured from 1D compression experiment

        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 462                                                       # Voxel units - vertical center
        xCenter = 507                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     		# Original GSD

    elif i == 4500 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-4500N/'
        ofl = '/home/eg/codes/pacOutput/1D/OTC-4500N/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details:
        dataName = 'otc-4500N'
        measVoidRatio = 0.359                                     # Void ratio measured from 1D compression experiment

        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 480                                                       # Voxel units - vertical center
        xCenter = 507                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     		# Original GSD


    eLen = 6*d50          # Edge length in mm

    # Reading and cropping the data file
    gliMap = Reader.readTiffFileSequence( inputFolderLocation,
                                          zCenter,
                                          yCenter,
                                          xCenter,
                                          eLen,
                                          cal,
                                          invImg=False)
    gsdOK = False

    # Naming tifffiles:
    gliName = ofl + 'gliMap.tiff'
    binName = ofl + 'binMap.tiff'
    edName = ofl + 'edMap.tiff'
    labName = ofl + 'labMap.tiff'
    corLabName = ofl + 'corLabMap.tiff'
    noEdgeCorLabName = ofl + 'noEdgeCorLabMap.tiff'

    while gsdOK == False:
        edmScaleUpFactor = int(input('Enter scaling for EDM: '))
        thresholdEdForPeak = int(input('Enter threshold of ED for peaks: '))

        binMap, binThresh, edMap, edPeakMap, labMap = Segment.obtLabMapITKWS( gliMap ,
                                                                              measuredVoidRatio=measVoidRatio ,
                                                                              outputLocation=ofl,
                                                                              edmScaleUp=edmScaleUpFactor,    # this represents how much the EDMs must be scaled up
                                                                              peakEdLimit=thresholdEdForPeak) # this represents what euclid distance should be considered for a peak

        
        corLabMap = Segment.fixErrSeg( labMap , pad=2, outputLocation=ofl , radiusRatioLimit=0.7)

        noEdgeCorLabMap = Segment.removeEdgeLabels( corLabMap )

        tf.imsave( gliName, gliMap.astype( 'uint32' ) )
        tf.imsave( binName, binMap.astype( 'uint32' ) )
        tf.imsave( edName, edMap.astype( 'uint32' ) )
        tf.imsave( labName, labMap.astype( 'uint32' ) )
        tf.imsave( corLabName , corLabMap.astype( 'uint32'))
        tf.imsave( noEdgeCorLabName , noEdgeCorLabMap.astype('uint32'))

        gsd1, gsd2, gsd3, gsd4, gsd5, gsd6= Measure.gsdAll( noEdgeCorLabMap , calib=cal )

        Plot.grainSizeDistribution(origGSD,gsd1,gsd2,gsd3,gsd4,gsd5,gsd6)

        exitLoop = input('\nIs any gsd ok(y/[n])?')
        #exitLoop = 'y'

        if exitLoop == 'y': gsdOK=True
        else: print('\n\nCheck the output file for (1) threshold error, (2) marker error')


    contactTableRW = Measure.contactNormalsSpam(corLabMap, method = 'rw')
    N, F, Fq = Measure.fabricVariablesWithUncertainity( contactTableRW, vectUncert = 0.26 )

    # Save files as csv
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd1.csv'), gsd1, delimiter=',')                        # Eqsp
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd2.csv'), gsd2, delimiter=',')                        # CA max
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd3.csv'), gsd3, delimiter=',')                        # CA med
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd4.csv'), gsd4, delimiter=',')                        # CA min
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd5.csv'), gsd5, delimiter=',')                        # Feret max
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd6.csv'), gsd6, delimiter=',')                        # Feret min
    np.savetxt((ofl+ str(eLen/d50) +'D50-contactTableRW.csv'), contactTableRW, delimiter=',')    # Contact table RW
    np.savetxt((ofl+ str(eLen/d50) +'N.txt'), N, fmt='%r')                                       # Fabric tensor
    np.savetxt((ofl+ str(eLen/d50) +'F.txt'), F, fmt='%r')                                       # Deviatoric fabric tensor
    np.savetxt((ofl+ str(eLen/d50) +'Fq.txt'), Fq, fmt='%r')                                     # Ansiotropy factor

    totalTimeEnd = time.time()
    totalTimeTaken = totalTimeEnd - totalTimeStart

    print('\n\n--------------------------------------**')
    print('Total time taken to analyze(mins): ~' + str(totalTimeTaken//60))