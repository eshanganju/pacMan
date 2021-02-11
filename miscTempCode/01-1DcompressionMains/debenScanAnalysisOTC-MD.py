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
from pac import Reader                      # Reads files, cuts into smaller pieces
from pac import Filter                      # Filters files using NLM filter
from pac import Segment                     # Binarizatio and WS
from pac import Measure                     # Calculates particle size, mprphology, contact, breakage
from pac import Plot                        # Plotting functions
import time                                 # The fourth dimension
import matplotlib.pyplot as plt             # matplotlib
import skimage.external.tifffile as tf      # scikit-image
import numpy as np                          # numpy
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
print('Starting')
totalTimeStart=time.time()

# 0 (0 MPa), 50 (1 MPa), 500 (10 MPa), 1500 (30 MPa)
data = np.array([0, 50, 500, 1500])

for i in data:
    if i == 0 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-0N/'
        ofl = '/home/eg/codes/pacOutput/OTC-MD-0N-7D50/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details 0N:
        dataName = 'otc-MD-0N'
        measuredVoidRatioSample = 0.609                                     # Void ratio measured from 1D compression experiment
        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 405                                                       # Voxel units - vertical center
        xCenter = 499                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

    elif i == 50 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-50N/'
        ofl = '/home/eg/codes/pacOutput/OTC-MD-50N-7D50/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details:
        dataName = 'otc-MD-50N'
        measuredVoidRatioSample = 0.60                                      # Void ratio measured from 1D compression experiment
        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 486                                                       # Voxel units - vertical center
        xCenter = 505                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

    elif i == 500 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-500N/'
        ofl = '/home/eg/codes/pacOutput/OTC-MD-500N-7D50/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details:
        dataName = 'otc-MD-500N'
        measuredVoidRatioSample = 0.585                                     # Void ratio measured from 1D compression experiment
        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 495                                                       # Voxel units - vertical center
        xCenter = 503                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

    elif i == 1500 :
        inputFolderLocation = '/home/eg/codes/pacInput/OTC-MD-1500N/'
        ofl = '/home/eg/codes/pacOutput/OTC-MD-1500N-7D50/'
        originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/otcOrig.csv' # Original GSD location

        # Data details:
        dataName = 'otc-MD-1500N'
        measuredVoidRatioSample = 0.56                                      # Void ratio measured from 1D compression experiment
        d50 = 0.72                                                          # D50 in mm - original gradation
        cal = 0.01193                                                       # calibration from CT mm/voxel
        zCenter = 513                                                       # Voxel units - center of slice
        yCenter = 497                                                       # Voxel units - vertical center
        xCenter = 505                                                       # Voxel units - horizontal center
        origGSD= np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

    eLen = 7*d50          # Edge length in mm

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
        # edmScaleUpFactor = int(input('Enter scaling for EDM: '))
        # thresholdEdForPeak = int(input('Enter threshold of ED for peaks: '))

        edmScaleUpFactor = 1
        thresholdEdForPeak = 5

        binMap, binThresh, edMap, edPeakMap, labMap = Segment.obtLabMapITKWS( gliMap ,
                                                                              measuredVoidRatio=measuredVoidRatioSample ,
                                                                              outputLocation=ofl,
                                                                              edmScaleUp = edmScaleUpFactor,           # this represents how much the EDMs must be scaled up
                                                                              peakEdLimit = thresholdEdForPeak)        # this represents what euclid distance should be considered for a peak

        corLabMap = Segment.fixErrSeg( labMap,
                                       pad=2,
                                       outputLocation=ofl,
                                       areaLimit=700 )

        '''
            Currently choosing areaLimit based on trial and error

            This can be some factor of the size of the particle and will depend on the resolution
            The diameters are 0.62, 0.72, 0.73 mm
            Average is 0.67mm
            That translates to 57 pixels
            Asuming a contact of half a particle i.e. 28 pixels, we get an area
            of 784 (square with edge 28 px)
            of 615 (circle with diameter 28 px)
            Average is around 700

            The area to be used should be a function of the sizes of the particles touching
            i.e. the contact between larger particles will be large and so for smaller
            This is especially true for crushed particles.
        '''

        noEdgeCorLabMap = Segment.removeEdgeLabels( corLabMap )
        gsd1, gsd2, gsd3, gsd4, gsd5, gsd6= Measure.gsdAll( noEdgeCorLabMap , calib=cal )

        #Plot.grainSizeDistribution(origGSD,gsd1,gsd2,gsd3,gsd4,gsd5,gsd6)

        tf.imsave( gliName, gliMap.astype( 'uint32' ) )
        tf.imsave( binName, binMap.astype( 'uint32' ) )
        tf.imsave( edName, edMap.astype( 'uint32' ) )
        tf.imsave( labName, labMap.astype( 'uint32' ) )
        tf.imsave( corLabName , corLabMap.astype( 'uint32'))
        tf.imsave( noEdgeCorLabName , noEdgeCorLabMap.astype('uint32'))

        #exitLoop = input('\nIs any gsd ok(y/[n])?')
        exitLoop = 'y'

        if exitLoop == 'y': gsdOK=True
        else: print('\n\nCheck the output file for (1) threshold error, (2) marker error')


    # Save files as csv
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd1.csv'), gsd1, delimiter=',')                        # Eqsp
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd2.csv'), gsd2, delimiter=',')                        # CA max
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd3.csv'), gsd3, delimiter=',')                        # CA med
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd4.csv'), gsd4, delimiter=',')                        # CA min
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd5.csv'), gsd5, delimiter=',')                        # Feret max
    np.savetxt((ofl+ str(eLen/d50) +'D50-gsd6.csv'), gsd6, delimiter=',')                        # Feret min

    # Contact table
    contactTableRW = Measure.contactNormalsSpam(corLabMap, method = 'rw')
    N, F, Fq = Measure.fabricVariablesWithUncertainity( contactTableRW, vectUncert = 0.26 )
    np.savetxt((ofl+ str(eLen/d50) +'D50-contactTableRW.csv'), contactTableRW, delimiter=',')    # Contact table RW
    np.savetxt((ofl+ str(eLen/d50) +'N.txt'), N, fmt='%r')    # Fabric tensor
    np.savetxt((ofl+ str(eLen/d50) +'F.txt'), F, fmt='%r')    # Deviatoric fabric tensor
    np.savetxt((ofl+ str(eLen/d50) +'Fq.txt'), Fq, fmt='%r')  # Ansiotropy factor

    totalTimeEnd = time.time()
    totalTimeTaken = totalTimeEnd - totalTimeStart

    print('\n\n--------------------------------------**')
    print('Total time taken to analyze(mins): ~' + str(totalTimeTaken//60))
