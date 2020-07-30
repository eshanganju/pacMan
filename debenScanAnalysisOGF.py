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
totalTimeStart=time.time()

# 0 (0 MPa), 100 (2 MPa), 500 (10 MPa), 1500 (30 MPa), 4500 (90 MPa)
data = 1500
maxPtclSize = 1         # mm
fracDim = 2.6     # Fractal dimension

if data == 0 :
    inputFolderLocation = '/home/eg/codes/pacInput/OGF-0N/'
    ofl = '/home/eg/codes/pacOutput/OGF-0N/' # output folder location ofl
    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

    # Data details 0N:
    dataName = 'ogf-0N'
    measuredVoidRatioSample = 0.635                                     # Void ratio measured from 1D compression experiment
    d50 = 0.62                                                          # D50 in mm - original gradation
    cal = 0.01193                                                       # calibration from CT mm/voxel
    zCenter = 508                                                       # Voxel units - center of slice
    yCenter = 437                                                       # Voxel units - vertical center
    xCenter = 490                                                       # Voxel units - horizontal center
    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

if data == 100 :
    inputFolderLocation = '/home/eg/codes/pacInput/OGF-100N/'
    ofl = '/home/eg/codes/pacOutput/OGF-100N/'
    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location


    # Data details 0N:
    dataName = 'ogf-100N'
    measuredVoidRatioSample = 0.625                                     # Void ratio measured from 1D compression experiment
    d50 = 0.62                                                          # D50 in mm - original gradation
    cal = 0.01193                                                       # calibration from CT mm/voxel
    zCenter = 508                                                       # Voxel units - center of slice
    yCenter = 447                                                       # Voxel units - vertical center
    xCenter = 490                                                       # Voxel units - horizontal center
    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

if data == 500 :
    inputFolderLocation = '/home/eg/codes/pacInput/OGF-500N/'
    ofl = '/home/eg/codes/pacOutput/OGF-500N/'
    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location


    # Data details 0N:
    dataName = 'ogf-500N'
    measuredVoidRatioSample = 0.615                                     # Void ratio measured from 1D compression experiment
    d50 = 0.62                                                          # D50 in mm - original gradation
    cal = 0.01193                                                       # calibration from CT mm/voxel
    zCenter = 508                                                       # Voxel units - center of slice
    yCenter = 458                                                       # Voxel units - vertical center
    xCenter = 490                                                       # Voxel units - horizontal center
    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

if data == 1500 :
    inputFolderLocation = '/home/eg/codes/pacInput/OGF-1500N/'
    ofl = '/home/eg/codes/pacOutput/OGF-1500N-9/'
    originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogf30MPa.csv' # Original GSD location


    # Data details 0N:
    dataName = 'ogf-1500N'
    measuredVoidRatioSample = 0.571  #0.565                             # Void ratio measured from 1D compression experiment
    d50 = 0.62                                                          # D50 in mm - original gradation
    cal = 0.01193                                                       # calibration from CT mm/voxel
    zCenter = 508                                                       # Voxel units - center of slice
    yCenter = 437                                                       # Voxel units - vertical center
    xCenter = 490                                                       # Voxel units - horizontal center
    origGSD = np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

eLen = 6*d50          # Edge length in mm

# Reading and cropping the data file
gliMap = Reader.readTiffFileSequence( inputFolderLocation, zCenter, yCenter, xCenter, eLen, cal, invertImageData=False)
gsdOK = False

# Save the 3D maps as tiff
gliName = ofl + 'gliMap.tiff'
binName = ofl + 'binMap.tiff'
edName = ofl + 'edMap.tiff'
labName = ofl + 'labMap.tiff'
correctedLabName = ofl + 'corLabMap.tiff'
noEdgeCorrectedLabName = ofl + 'noEdgeCorLabMap.tiff'

while gsdOK == False:
    binMap, edMap, edPeakMap, labMap = Segment.obtainLabelledMapUsingITKWS( gliMap , measuredVoidRatio=measuredVoidRatioSample , outputLocation=ofl )

    correctedLabMap = Segment.fixErrorsInSegmentation( labMap , pad=2, outputLocation=ofl , areaLimit = 650)
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
    #correctedLabMap = labMap

    noEdgeCorrectedLabMap = Segment.removeEdgeLabels( correctedLabMap )

    gsd1, gsd2, gsd3, gsd4 = Measure.gsd( noEdgeCorrectedLabMap , calib=cal )

    #junk,junk,junk,Br1 = Measure.relBreak(origGSD,gsd1, maxSize=maxPtclSize, fracDim=fracDim)
    #junk,junk,junk,Br2 = Measure.relBreak(origGSD,gsd2, maxSize=maxPtclSize, fracDim=fracDim)
    #junk,junk,junk,Br3 = Measure.relBreak(origGSD,gsd3, maxSize=maxPtclSize, fracDim=fracDim)
    #junk,junk,junk,Br4 = Measure.relBreak(origGSD,gsd4, maxSize=maxPtclSize, fracDim=fracDim)

    #print('Br1 = ' + str(Br1))
    #print('Br2 = ' + str(Br2))
    #print('Br3 = ' + str(Br3))
    #print('Br4 = ' + str(Br4))

    Plot.grainSizeDistribution(origGSD,gsd1,gsd2,gsd3,gsd4)

    tf.imsave( gliName, gliMap.astype( 'uint32' ) )
    tf.imsave( binName, binMap.astype( 'uint32' ) )
    tf.imsave( edName, edMap.astype( 'uint32' ) )
    tf.imsave( labName, labMap.astype( 'uint32' ) )
    tf.imsave( correctedLabName , correctedLabMap.astype( 'uint32'))
    tf.imsave( noEdgeCorrectedLabName , noEdgeCorrectedLabMap.astype('uint32'))

    exitLoop = input('\nIs any gsd ok(y/[n])?')

    if exitLoop == 'y':
        gsdOK=True
        #gsdNum = int(input('Which gsd is best for Br calcs(1,2,3,4)?:'))
        gsdNum = 4

        if gsdNum == 1 : gsdForBr = gsd1
        elif gsdNum == 2 : gsdForBr = gsd2
        elif gsdNum == 3 : gsdForBr = gsd3
        elif gsdNum == 4 : gsdForBr = gsd4

    else: print('\nUse user threshold to update binary map - check the binaryThreshold file in output folder for latest threshold')

formatGsdOrig, formatGsdCurr, formatGsdUlt, Br = Measure.relBreak(origGSD,gsdForBr)

contactTableRW = Measure.contactNormalsSpam(noEdgeCorrectedLabMap, method = 'rw')
N, F, Fq = Measure.fabricVariablesWithUncertainity( contactTableRW, vectUncert = 0.26 )
Plot.equalAreaProjection(contactTableRW)

# Save files as csv
np.savetxt((ofl+ str(eLen/d50) +'D50-gsd1.csv'), gsd1, delimiter=',')                        # Eqsp
np.savetxt((ofl+ str(eLen/d50) +'D50-gsd2.csv'), gsd2, delimiter=',')                        # CA max
np.savetxt((ofl+ str(eLen/d50) +'D50-gsd3.csv'), gsd3, delimiter=',')                        # CA med
np.savetxt((ofl+ str(eLen/d50) +'D50-gsd4.csv'), gsd4, delimiter=',')                        # CA min
np.savetxt((ofl+ str(eLen/d50) +'D50-contactTableRW.csv'), contactTableRW, delimiter=',')    # Contact table RW
np.savetxt((ofl+ str(eLen/d50) +'N.txt'), N, fmt='%r')    # Contact table RW
np.savetxt((ofl+ str(eLen/d50) +'F.txt'), F, fmt='%r')    # Contact table RW
np.savetxt((ofl+ str(eLen/d50) +'Fq.txt'), Fq, fmt='%r')    # Contact table RW

brFile = open(ofl+ str(eLen) +'D50-Br.txt',"w")
L = 'Br = ' + str(Br) + '%'
brFile.write(L)
brFile.close()

totalTimeEnd = time.time()
totalTimeTaken = totalTimeEnd - totalTimeStart

print('\n\n--------------------------------------**')
print('Total time taken to analyze(mins): ~' + str(totalTimeTaken//60))
