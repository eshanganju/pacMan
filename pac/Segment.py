'''
Description:
    Segment module of PAC.
'''

from scipy.ndimage.morphology import distance_transform_edt as edt
from scipy.ndimage.morphology import binary_fill_holes
from scipy.ndimage.morphology import binary_opening
from skimage.morphology import local_maxima as localMaxima
from skimage.morphology import h_maxima as hmax
from skimage.morphology import watershed as wsd
from skimage.filters import threshold_otsu
import tifffile as tiffy
import matplotlib.pyplot as plt
import numpy as np
import spam.label as slab
import time
import math
from numba import jit

from pac import Measure

# This is to plot all the text when running the functions
VERBOSE = True
TESTING = True

def _obtLabMapITKWS( gliMap , knownThreshold = None, measuredVoidRatio = None, outputLoc=None, edmScaleUp=1, peakEdLimit=5):
    """
    Description:

    Parameters:

    Return:
    """
    print( '\nSegmenting particles by ITK topological watershed' )

    if knownThreshold == None:
        if measuredVoidRatio == None:
            binMethod=input('Which binarization method to use (1) OTSU; (2) User; [3] Density?: ')
        else:
            binMethod = '3'
        if binMethod == '1': binThresh, binMask = binarizeAccordingToOtsu( gliMap )
        elif binMethod == '2': binThresh, binMask = binarizeAccordingToUserThreshold( gliMap  )
        else : binThresh, binMask = binarizeAccordingToDensity( gliMap , measuredVoidRatio )

    else :
        print('User threshold with :' + str(round(knownThreshold)))
        binThresh, binMask = binarizeAccordingToUserThreshold( gliMap, knownThreshold )

    voidRatio = calcVoidRatio( binMask )

    if outputLoc != None:
        threshFile = open(outputLoc + 'binaryThreshold.txt',"a+")
        threshFile.write( "Bin Threshold = " + str( binThresh ) )
        threshFile.write( "\nVoid ratio = " + str( voidRatio ) + '\n\n')
        threshFile.close()

    # Simple Euclidean distance
    edMap = obtainEuclidDistanceMap( binMask , scaleUp=edmScaleUp)

    # Peaks in EDM
    edPeaksMap = obtainLocalMaximaMarkers( edMap , h=peakEdLimit)

    print('\n\nStarting ITK WS')
    print( '--------------------------*' )
    labelledMap = wsd( -edMap, markers = edPeaksMap, mask = binMask )
    print( 'Watershed segmentation complete' )

    return binMask, binThresh, edMap, edPeaksMap, labelledMap

def segmentUsingWatershed(binaryMapToSeg,edmMapForTopo,edmPeaksForSeed,sampleName='',saveImg=True,outputDir=''):
    """Simple function that uses skimage watershed and saves a copy of the segmented image

    Parameters
    ----------
    binaryMapToSeg : ndarray

    edmMapForTopo : ndarray

    edmPeaksForSeed : ndarray

    sampleName : string

    saveImg : bool

    outputDir : string

    Return
    ------
    labMap : ndarray
        labelled map with each particle assigned a separate integer value
    """
    print('\nStarting segmentation using watershed')
    print('-----------------------------------------*')
    labMap = wsd(-edmMapForTopo,markers=edmPeaksForSeed,mask=binaryMapToSeg)

    if saveImg == True:
        if VERBOSE: print('\nSaving lab map...')
        tiffy.imsave(outputDir+sampleName+'-labMap.tif',labMap.astype('uint16'))

    return labMap.astype('uint16')

def binarizeAccordingToOtsu( gliMapToBinarize, sampleName='', saveImg=False, outputDir='', returnThresholdVal=False):
    """
    Description:
        Function to binarize GLI map according to OTSUs algorithm
        Uses the skimage.filter module threshold_otsu

    Parameters:
        gliMapToBinarize (nd array): (preferably) filtered GLI Map
        fileName (str): name of the sample
        saveImg (bool): Should we save the image or not
        outputLic (str): Location of the output image

    Returns:
        otsuThreshold (float): duh
        binaryMap (nd array): Array binarized.

    """
    print('\nRunning Otsu Binarization')
    print('----------------------------*')
    otsuThreshold = threshold_otsu( gliMapToBinarize )
    binaryMap = np.zeros_like( gliMapToBinarize )

    binaryMap[ np.where( gliMapToBinarize > otsuThreshold ) ] = 1
    binaryMap = binaryMap.astype( int )
    e1 = calcVoidRatio( binaryMap )

    binaryMap = fillHoles( binaryMap )
    e2 = calcVoidRatio( binaryMap )

    binaryMap = removeSpecks( binaryMap )
    e3 = calcVoidRatio( binaryMap )

    if VERBOSE == True: print( 'Global Otsu threshold = %f' % otsuThreshold )
    if VERBOSE == True: print( 'Void ratio after threshold = %f' % e1 )
    if VERBOSE == True: print( 'Void ratio after filling holes = %f' % e2 )
    if VERBOSE == True: print( 'Void ratio after removing specks = %f' % e3 )

    if saveImg == True:
        if VERBOSE: print('\nSaving binary map...')
        tiffy.imsave( outputDir + sampleName + '-binaryMap.tif', binaryMap.astype('uint16'))
    otsuThresholdFileName = outputDir+sampleName+'-otsuThreshold.txt'

    f = open( otsuThresholdFileName,"w+" )
    f.write( 'Otsus threshold = %f' % otsuThreshold )
    f.close()

    if returnThresholdVal == True: return otsuThreshold, binaryMap

    elif returnThresholdVal == False: return binaryMap

def binarizeAccordingToUserThreshold( gliMapToBinarize, userThreshold = None, returnThresholdVal=False ):
    """
    Description:
        Binarize according to user supplied threshold

    Parameters:
        gliMapToBinarize
        userThreshold

    Return:
        Binarized map and user threshold
    """
    if userThreshold == None:
        userThreshold = int( input( 'Enter user threshold: ' ) )

    binaryMap = np.zeros_like( gliMapToBinarize )
    binaryMap[ np.where( gliMapToBinarize > userThreshold ) ] = 1
    binaryMap = binaryMap.astype( int )
    e1 = calcVoidRatio( binaryMap )

    # Filling Holes
    binaryMap = fillHoles( binaryMap )
    e2 = calcVoidRatio( binaryMap )

    # Removing specks
    binaryMap = removeSpecks( binaryMap )
    e3 = calcVoidRatio( binaryMap )

    print( '\nGlobal User threshold = %f' % userThreshold )
    print( 'Void ratio after threshold = %f' % e1 )
    print( 'Void ratio after filling holes = %f' % e2 )
    print( 'Void ratio after removing specks = %f' % e3 )

    if returnThresholdVal == True: return userThreshold, binaryMap
    if returnThresholdVal == False: return binaryMap

def binarizeAccordingToDensity( gliMapToBinarize , measuredVoidRatio = None):
    print('\nRunning density-based threshold...')
    print('Running Otsu first to get inital guess: ')
    otsuThreshold, otsuBinMap = binarizeAccordingToOtsu( gliMapToBinarize )
    currentThreshold = otsuThreshold

    if measuredVoidRatio == None: measuredVoidRatio = int( input('Input the known void ratio: ') )

    currentVoidRatio = calcVoidRatio( otsuBinMap )
    targetVoidRatio = measuredVoidRatio
    greyLvlMap = gliMapToBinarize

    tolerance = 0.001
    deltaVoidRatio = targetVoidRatio - currentVoidRatio

    absDeltaThreshold = 0
    maxAbsDeltaThreshold = 500

    iterationNum = 1
    maxIterations = 50

    print('\nRunning iterations to compute threshold corresponding to target void ratio...')
    print('Target void ratio = ' + str( targetVoidRatio ) )

    while( abs( deltaVoidRatio ) > tolerance and iterationNum <= maxIterations ):
        incrementSign = deltaVoidRatio/abs(deltaVoidRatio)  

        if abs( deltaVoidRatio ) > 0.100:
            absDeltaThreshold = 300

        elif abs( deltaVoidRatio ) > 0.05:
            absDeltaThreshold = 150

        elif abs( deltaVoidRatio ) > 0.01:
            absDeltaThreshold = 50

        elif abs( deltaVoidRatio ) > 0.005:
            absDeltaThreshold = 25

        elif abs( deltaVoidRatio ) > 0.001:
            absDeltaThreshold = 10

        elif abs( deltaVoidRatio ) > 0.0005:
            absDeltaThreshold = 5

        else:
            absDeltaThreshold = 1

        currentThreshold = currentThreshold + absDeltaThreshold*incrementSign
        currentBinaryMap = np.zeros_like( greyLvlMap )        
        currentBinaryMap[ np.where( greyLvlMap > currentThreshold ) ] = 1

        currentBinaryMap = fillHoles( currentBinaryMap )
        currentBinaryMap = removeSpecks( currentBinaryMap )

        currentVoidRatio = calcVoidRatio( currentBinaryMap )
        deltaVoidRatio = targetVoidRatio - currentVoidRatio

        print( '\nIteration %d:' % iterationNum )
        print( 'Current threshold: %d' % round( currentThreshold ) )
        print( 'Otsu Threshold: %d' % round( otsuThreshold ) )
        print( 'Current void ratio: %0.5f' % round( currentVoidRatio,5 ) )
        print( 'Target void ratio: %0.5f' % round( targetVoidRatio,5 ) )
        print( 'Target - Current void ratio: %0.5f' % round( deltaVoidRatio, 5 ) )

        #f = open( userThresholdStepsTextFileName, 'a' )
        #f.write( "\nIteration %d:\n" % iterationNum)
        #f.write( "Current threshold: %d\n" % round( currentThreshold ) )
        #f.write( "Otsu Threshold: %d\n" %  round( otsuThreshold ) )
        #f.write( "Current void ratio:  %0.5f\n" % round( currentVoidRatio,5 ) )
        #f.write( "Target void ratio: %0.5f\n" % round( targetVoidRatio,5 ) )
        #f.write( "Target - Current void ratio: %0.5f\n" % round( deltaVoidRatio, 5 ) )
        #f.close()

        iterationNum = iterationNum + 1

    #densityBasedThresholdTextFileName = outputLoc +  sampleName + '-userThresholdDensityBased.txt'
    #f = open( densityBasedThresholdTextFileName, "w+" )
    #f.write( "Global density based threshold = %f\n" % currentThreshold )
    #f.close()

    print('\nThreshold corresponding to measured density = ' + str( np.round( currentThreshold ) ) )
    print('-------------------------------------------------------------*')

    return currentThreshold, currentBinaryMap

def calcVoidRatio( binaryMapforVoidRatioCalc ):
    '''
    Description:
        Caculated void ratio from binary map

    Parameters:
        binaryMap (int array): numpy array with 1 and 0

    Returns
        void ratio of the binary map assumin 1 is particle
    '''
    volSolids = binaryMapforVoidRatioCalc.sum()
    lenZ = binaryMapforVoidRatioCalc.shape[ 0 ]
    lenY = binaryMapforVoidRatioCalc.shape[ 1 ]
    lenX = binaryMapforVoidRatioCalc.shape[ 2 ]
    volTotal =  lenZ * lenY * lenX 
    currentVoidRatio = ( volTotal - volSolids ) / volSolids
    return currentVoidRatio

def fillHoles( oldBinaryMapWithHoles ):
    allOk = False
    newNoHoleBinaryMap=oldBinaryMapWithHoles
    newNoHoleBinaryMap = binary_fill_holes( newNoHoleBinaryMap )
    return newNoHoleBinaryMap.astype(int)

def removeSpecks( oldBinaryMapWithSpecks ):
    remSpecksCheck = 'y'

    if remSpecksCheck == 'y':newNoSpekBinaryMap = binary_opening( oldBinaryMapWithSpecks )
    else: newNoSpekBinaryMap = oldBinaryMapWithSpecks

    return newNoSpekBinaryMap.astype(int)

def obtainEuclidDistanceMap( binaryMapForEDM, scaleUp = int(1), saveImg=False, sampleName='', outputDir=''):
    """Computes the euclidian distance tranform (EDT) for a binary map

    An EDT or an euclidian distance map (EDM) is an ndarray of the same size as
    the input binary map. Instead of 1 or 0, as in a binary map, each solid (1)
    voxel is assigned the distance between it and the nearest void (0) voxel.
    The resulting map shows how far each voxel is to the nearest void voxel,
    i.e. from the boundary of paricle.

    Parameters
    ----------
    binaryMapForEDM : ndarray
        binary map of the xct scan

    scaleUp : unsigned integer
        This is done inorder to artificially increase the zoom of the data. The
        distance from the voxels increased if a fractional value is needed in
        the location of peaks.

    saveImg : bool
        If this is true, the edt is saved in the requested location or at the
        location whenere the code is run

    sampleName : string
        name of the sample - used to name the file used to store the edt data

    outputDir : string
        Location of the the ouput directory. If left empty, the file is saved at
        the same location as the location of the code.

    Returns
    -------
    edMap : ndarray
        array containing the euclidian distance map of the binary map

    """
    print('\nFinding Euclidian distance map (EDM)')
    print('------------------------------------*')

    edMap = edt( binaryMapForEDM )
    if scaleUp!=0 : edMap =  edMap * scaleUp

    print( "EDM Created" )

    if saveImg == True:
        if VERBOSE: print('\nSaving EDM map...')
        tiffy.imsave(outputDir + sampleName + '-edm.tif',edMap)

    return edMap

def obtainLocalMaximaMarkers( edMapForPeaks , method = 'hlocal' , h=5, saveImg=False, sampleName='', outputDir=''):
    """Computes the local maximas in the euclidean distance map.
    it uses skimage.morphology

    Parameters
    ----------
       edMapForPeaks : ndarray containing the euclidian distance map.
       method : string containing the choice of algorithm
       h=5,
       saveImg=False,
       sampleName='',
       outputDir=''

    Returns
    -------
    edmPeakMarkers
        ndArray of the same size as the input array containing map
        of peaks in the edm, numbered in an increasing order
    """
    if VERBOSE == True:
        print('\nObtaining peaks of EDM...')
        print('------------------------------*')
        print( 'Finding local maxima in EDM' )

    if method == 'hlocal' :
        if h == None: h = int( input( 'Enter the minimum height for a peak (px): ') )
        edmPeakMarkers = hmax( edMapForPeaks, h ).astype(int)

    elif method == 'local' : edmPeakMarkers = localMaxima( edMapForPeaks ).astype(int)

    if VERBOSE:
        print( '\tFound local maximas' )
        print( '\nResetting count of peaks' )

    # Resetting counts such that each peak has a unique integer label
    count=0
    for frame in range( 0, edmPeakMarkers.shape[ 0 ] ):
        for row in range( 0, edmPeakMarkers.shape[ 1 ] ):
            for col in range( 0, edmPeakMarkers.shape[ 2 ] ):
                if edmPeakMarkers[ frame ][ row ][ col ] == 1:
                    edmPeakMarkers[ frame ][ row ][ col ] = count + 1
                    count = count + 1
        if VERBOSE:
            print( 'Processed ' + str(frame + 1) + ' out of ' + str(edmPeakMarkers.shape[ 0 ]) + ' slices' )

    if VERBOSE:
        print('\tCounts reset')
        print('\tNumber of peaks: ' + str(round(edmPeakMarkers.max())))

    if saveImg == True:
        if VERBOSE: print('\nSaving EDM peaks map...')
        tiffy.imsave( outputDir + sampleName + '-edmPeaks.tif', edmPeakMarkers.astype('uint16'))

    return edmPeakMarkers

@jit(nopython=True)
def fixErrorsInSegmentationWithNumba( labelledMapForOSCorr, pad=2, areaLimit = 700,
                                      conside rEdgeLabels=True,checkForSmallParticles = True,
                                      radiusCheck=True, radiusRatioLimit=0.5, sampleName='',
                                      saveImg=True, outputDir=''):
    """Corrects over segmentation caused by incorrect edm peak selection

    There are two main approaches. One way is to compute the "area of contact"
    between the two particles suspected of being one, and if the area is larger
    than an absolute threshold, then the two particles are merged

    The other approach is to compute the contact area and convert it to an
    equivalent radii (assuming a circular contact area). If the ratio of this
    radii (of the contact) to the radii of the particles (assuming an equivanent
    sphere) is greater than a threshold, then the two particles are merged

    The contact are a is determined using spam.contact's contactingLabels
    function.

    Parameters:
        labelledMapForOSCorr : ndarray
        pad : unsigned integer
        areaLimit : unsigned integer
        considerEdgeLabels : bool
        checkForSmallParticles : bool
        radiusCheck : bool
        radiusRatioLimit : float
        sampleName : string
        saveImg : bool
        outputDir : string

    Return:
        correctedCleanedLabelMap : ndarray
            Corrected label map
    """

    print('\nStarting label correction')
    print('---------------------------*')
    #if areaLimit != None : print('Area limit is : ' + str( np.round( areaLimit) ) )

    if radiusCheck == True: print('Radius ratio limit is: ' + str( radiusRatioLimit) )
    elif radiusCheck == False: print('Area limit is: ' + str(areaLimit))

    # Apply padding to the data
    if pad > 0:
        labelledMap = labelledMapForOSCorr
        padLabMap = np.zeros( ( labelledMap.shape[0]+2*pad, labelledMap.shape[0]+2*pad, labelledMap.shape[0]+2*pad ) )
        padLabMap[pad : padLabMap.shape[0]-pad , pad : padLabMap.shape[1]-pad , pad : padLabMap.shape[ 2 ]-pad ] = labelledMap
        labelledMapForOSCorr = padLabMap

    lastLabel = labelledMapForOSCorr.max()
    currentLabel = 1
    correctedLabelMap = labelledMapForOSCorr
    timeStart = time.time()

    if VERBOSE: print('Currently padding is ' + str(pad) + ' px')

    # Include the edge labels in the calculation
    if considerEdgeLabels == False:
        if VERBOSE: print('\tOk, not consolidating edge labels\n')
        edgePad = pad
        # The edgepad adds additional padding to prevent the edge lables from being removed
    else:
        if VERBOSE: print('\tOk, consolidating edge labels\n')
        edgePad = 0

    if outputDir != '' : correctionLog = open(outputDir + sampleName + '-correctionLog.Txt',"a+")
    else: correctionLog = open("correctionLog.txt","a+")

    if VERBOSE:
        print('Starting edge label consolidation - This may take some time')

    correctionLog.write('\n-------------------------------------------------------')

    if radiusCheck == False: correctionLog.write('\nArea limit: ' + str( np.round( areaLimit ) ) + '\n\n')
    elif radiusCheck == True: correctionLog.write('\nRadius ratio limit: ' + str( radiusRatioLimit ) + '\n\n')

    # Loop through labels till all OS labels are merged
    while currentLabel <= lastLabel:
        isEdgeLabel = checkIfEdgeLabel( correctedLabelMap, currentLabel, edgePad )

        # If edge label, move to next label
        if isEdgeLabel == True:
            if VERBOSE: print('Label %d is an edge label, moving to next label' % currentLabel)
            correctionLog.write('\n\nLabel ' + str(currentLabel) + ' is an edge label, moving to next label')
            currentLabel += 1

        # If not edge label, check contacting labels
        else:
            contactLabel, contactArea = slab.contactingLabels(correctedLabelMap, currentLabel, areas=True)
            if VERBOSE: print('\n')

            currentLabelRadius = Measure.getEqspDia(correctedLabelMap,currentLabel)[1]/2
            touchingParticleRadius = []
            contactRadius = []
            radiusRatio =[]

            for touchingLabel in contactLabel:
                touchingParticleRadius.append((Measure.getEqspDia(correctedLabelMap,touchingLabel)[1])/2 )

            for touchingArea in contactArea:
                contactAreaRadiusVal = 0.5 * ( 4 * touchingArea / math.pi ) ** 0.5
                contactRadius.append(contactAreaRadiusVal)

            for contact in range( 0 , len(contactLabel) ):
                minR = min( currentLabelRadius,touchingParticleRadius[contact])
                r = contactRadius[contact]
                radiusRatio.append(r/minR)

            largeAreaVal = [contactArea[location] for location, val in enumerate(contactArea) if val >= areaLimit]
            largeAreaLabel = [contactLabel[location] for location, val in enumerate(contactArea) if val >= areaLimit]

            largeRatioVal = [radiusRatio[location] for location, ratio in enumerate(radiusRatio) if ratio >= radiusRatioLimit]
            largeRatioLabel = [contactLabel[location] for location, ratio in enumerate(radiusRatio) if ratio >= radiusRatioLimit]

            if VERBOSE: print('\nLabel ' + str( currentLabel ) + ' is contacting ' + str( contactLabel ) )
            if VERBOSE: print('With radius Ratios: ' + str( radiusRatio ) )

            correctionLog.write('\n\nLabel ' + str( currentLabel ) + ' is contacting ' + str( contactLabel ) )
            correctionLog.write('\nWith areas: ' + str( contactArea ) )
            correctionLog.write('\nCurrent particle radius: ' + str( currentLabelRadius ) )
            correctionLog.write('\nTouching particles radius: ' + str( touchingParticleRadius ) )
            correctionLog.write('\nRadius ratios: ' + str( radiusRatio ) )

            # Merging with radius ratio limits
            if radiusCheck == True:
                if len(largeRatioVal) != 0:
                    if VERBOSE: print('\tLabels with large ratios: ' + str( largeRatioLabel ) )
                    if VERBOSE: print('\tRatio for above labels: ' + str( largeRatioVal ) )

                    correctionLog.write('\n\tLabels with large ratio: ' + str( largeRatioLabel ) )
                    correctionLog.write('\n\tRatio for above labels: ' + str( largeRatioVal ) )

                    labelToMerge = largeRatioLabel[ largeRatioVal.index( min( largeRatioVal ) ) ]
                    smallerLabel = min( labelToMerge, currentLabel )
                    largerLabel = max( labelToMerge, currentLabel )
                    correctedLabelMap[np.where(correctedLabelMap == largerLabel)] = int(smallerLabel)
                    correctedLabelMap = moveLabelsUp( correctedLabelMap,largerLabel )
                    currentLabel = smallerLabel
                    lastLabel = correctedLabelMap.max()

                    if VERBOSE: print('\tMerging label %d and %d' %(smallerLabel, largerLabel))
                    correctionLog.write('\n\tMerging labels' + str(smallerLabel) + ' and ' + str(largerLabel))
                    if VERBOSE: print( 'Checking from label %d again.' % currentLabel )
                    correctionLog.write('\nChecking from label ' + str(currentLabel ) + ' again.')

                else:
                    if VERBOSE: print('Label ' + str(currentLabel) + ' has no large contacts.')
                    correctionLog.write('\nLabel ' + str(currentLabel) + ' has no large contacts.')
                    currentLabel = currentLabel + 1
                    if VERBOSE: print( 'Moving to next label ' + str(currentLabel) +  ' now.')
                    correctionLog.write('\nMoving to next label ' + str(currentLabel) +  ' now.')

            # Merging with area limits
            else:
                if len(largeAreaVal) != 0:
                    if VERBOSE: print('\tLabels with large area: ' + str( largeAreaLabel ) )
                    if VERBOSE: print('\tArea for above labels: ' + str( largeAreaVal ) )

                    correctionLog.write('\n\tLabels with large area: ' + str( largeAreaLabel ) )
                    correctionLog.write('\n\tArea for above labels: ' + str( largeAreaVal ) )

                    labelToMerge = largeAreaLabel[ largeAreaVal.index( min( largeAreaVal ) ) ]
                    smallerLabel = min(labelToMerge,currentLabel)
                    largerLabel = max(labelToMerge,currentLabel)
                    correctedLabelMap[np.where(correctedLabelMap == largerLabel)] = int(smallerLabel)
                    correctedLabelMap = moveLabelsUp( correctedLabelMap , largerLabel )
                    currentLabel = smallerLabel
                    lastLabel = correctedLabelMap.max()

                    if VERBOSE: print('\tMerging label %d and %d' %(smallerLabel, largerLabel))
                    correctionLog.write('\n\tMerging labels ' + str(smallerLabel) + ' and ' + str(largerLabel))
                    if VERBOSE: print( 'Checking from label %d again.' % currentLabel )
                    correctionLog.write('\nChecking from label ' + str(currentLabel ) + ' again.')

                else:
                    if VERBOSE: print('Label ' + str(currentLabel) + ' has no large contacts.')
                    correctionLog.write('\n\nLabel ' + str(currentLabel) + ' has no large contacts.')
                    currentLabel = currentLabel + 1
                    if VERBOSE: print( 'Moving to next label ' + str(currentLabel) +  ' now.')
                    correctionLog.write('\nMoving to next label ' + str(currentLabel) +  ' now.')

    if pad > 0: correctedLabelMap = removePaddingFromLabelledMap(correctedLabelMap, pad)

    timeEnd = time.time()
    timeTaken = (timeEnd - timeStart)//60

    print( '\nTime taken for correction loop: ' + str( np.round(timeTaken) ) + ' mins' )
    correctionLog.write('\nTime taken for correction loop: ' + str( np.round(timeTaken) ) + ' mins\n\n')

    correctionLog.close()

    if checkForSmallParticles == True:
        correctedCleanedLabelMap = removeSmallParticles( correctedLabelMap )
    else : correctedCleanedLabelMap = correctedLabelMap

    if saveImg == True:
        print('\nSaving corrected labelled map...')
        tiffy.imsave(outputDir+sampleName+'-correctedLabelMap.tif',correctedCleanedLabelMap.astype('uint16'))

    return correctedCleanedLabelMap

def fixErrorsInSegmentation( labelledMapForOSCorr, pad=2, areaLimit = 700,
                             conside rEdgeLabels=True,checkForSmallParticles = True,
                             radiusCheck=True, radiusRatioLimit=0.5, sampleName='',
                             saveImg=True, outputDir=''):
    """Corrects over segmentation caused by incorrect edm peak selection

    There are two main approaches. One way is to compute the "area of contact"
    between the two particles suspected of being one, and if the area is larger
    than an absolute threshold, then the two particles are merged

    The other approach is to compute the contact area and convert it to an
    equivalent radii (assuming a circular contact area). If the ratio of this
    radii (of the contact) to the radii of the particles (assuming an equivanent
    sphere) is greater than a threshold, then the two particles are merged

    The contact are a is determined using spam.contact's contactingLabels
    function.

    Parameters:
        labelledMapForOSCorr : ndarray
        pad : unsigned integer
        areaLimit : unsigned integer
        considerEdgeLabels : bool
        checkForSmallParticles : bool
        radiusCheck : bool
        radiusRatioLimit : float
        sampleName : string
        saveImg : bool
        outputDir : string

    Return:
        correctedCleanedLabelMap : ndarray
            Corrected label map
    """

    print('\nStarting label correction')
    print('---------------------------*')
    #if areaLimit != None : print('Area limit is : ' + str( np.round( areaLimit) ) )

    if radiusCheck == True: print('Radius ratio limit is: ' + str( radiusRatioLimit) )
    elif radiusCheck == False: print('Area limit is: ' + str(areaLimit))

    # Apply padding to the data
    if pad > 0: labelledMapForOSCorr = applyPaddingToLabelledMap(labelledMapForOSCorr, pad)

    lastLabel = labelledMapForOSCorr.max()
    currentLabel = 1
    correctedLabelMap = labelledMapForOSCorr
    timeStart = time.time()

    if VERBOSE: print('Currently padding is ' + str(pad) + ' px')

    # Include the edge labels in the calculation
    if considerEdgeLabels == False:
        if VERBOSE: print('\tOk, not consolidating edge labels\n')
        edgePad = pad
        # The edgepad adds additional padding to prevent the edge lables from being removed
    else:
        if VERBOSE: print('\tOk, consolidating edge labels\n')
        edgePad = 0

    if outputDir != '' : correctionLog = open(outputDir + sampleName + '-correctionLog.Txt',"a+")
    else: correctionLog = open("correctionLog.txt","a+")

    if VERBOSE:
        print('Starting edge label consolidation - This may take some time')

    correctionLog.write('\n-------------------------------------------------------')

    if radiusCheck == False: correctionLog.write('\nArea limit: ' + str( np.round( areaLimit ) ) + '\n\n')
    elif radiusCheck == True: correctionLog.write('\nRadius ratio limit: ' + str( radiusRatioLimit ) + '\n\n')

    # Loop through labels till all OS labels are merged
    while currentLabel <= lastLabel:
        isEdgeLabel = checkIfEdgeLabel( correctedLabelMap, currentLabel, edgePad )

        # If edge label, move to next label
        if isEdgeLabel == True:
            if VERBOSE: print('Label %d is an edge label, moving to next label' % currentLabel)
            correctionLog.write('\n\nLabel ' + str(currentLabel) + ' is an edge label, moving to next label')
            currentLabel += 1

        # If not edge label, check contacting labels
        else:
            contactLabel, contactArea = slab.contactingLabels(correctedLabelMap, currentLabel, areas=True)
            if VERBOSE: print('\n')

            currentLabelRadius = Measure.getEqspDia(correctedLabelMap,currentLabel)[1]/2
            touchingParticleRadius = []
            contactRadius = []
            radiusRatio =[]

            for touchingLabel in contactLabel:
                touchingParticleRadius.append((Measure.getEqspDia(correctedLabelMap,touchingLabel)[1])/2 )

            for touchingArea in contactArea:
                contactAreaRadiusVal = 0.5 * ( 4 * touchingArea / math.pi ) ** 0.5
                contactRadius.append(contactAreaRadiusVal)

            for contact in range( 0 , len(contactLabel) ):
                minR = min( currentLabelRadius,touchingParticleRadius[contact])
                r = contactRadius[contact]
                radiusRatio.append(r/minR)

            largeAreaVal = [contactArea[location] for location, val in enumerate(contactArea) if val >= areaLimit]
            largeAreaLabel = [contactLabel[location] for location, val in enumerate(contactArea) if val >= areaLimit]

            largeRatioVal = [radiusRatio[location] for location, ratio in enumerate(radiusRatio) if ratio >= radiusRatioLimit]
            largeRatioLabel = [contactLabel[location] for location, ratio in enumerate(radiusRatio) if ratio >= radiusRatioLimit]

            if VERBOSE: print('\nLabel ' + str( currentLabel ) + ' is contacting ' + str( contactLabel ) )
            if VERBOSE: print('With radius Ratios: ' + str( radiusRatio ) )

            correctionLog.write('\n\nLabel ' + str( currentLabel ) + ' is contacting ' + str( contactLabel ) )
            correctionLog.write('\nWith areas: ' + str( contactArea ) )
            correctionLog.write('\nCurrent particle radius: ' + str( currentLabelRadius ) )
            correctionLog.write('\nTouching particles radius: ' + str( touchingParticleRadius ) )
            correctionLog.write('\nRadius ratios: ' + str( radiusRatio ) )

            # Merging with radius ratio limits
            if radiusCheck == True:
                if len(largeRatioVal) != 0:
                    if VERBOSE: print('\tLabels with large ratios: ' + str( largeRatioLabel ) )
                    if VERBOSE: print('\tRatio for above labels: ' + str( largeRatioVal ) )

                    correctionLog.write('\n\tLabels with large ratio: ' + str( largeRatioLabel ) )
                    correctionLog.write('\n\tRatio for above labels: ' + str( largeRatioVal ) )

                    labelToMerge = largeRatioLabel[ largeRatioVal.index( min( largeRatioVal ) ) ]
                    smallerLabel = min( labelToMerge, currentLabel )
                    largerLabel = max( labelToMerge, currentLabel )
                    correctedLabelMap[np.where(correctedLabelMap == largerLabel)] = int(smallerLabel)
                    correctedLabelMap = moveLabelsUp( correctedLabelMap,largerLabel )
                    currentLabel = smallerLabel
                    lastLabel = correctedLabelMap.max()

                    if VERBOSE: print('\tMerging label %d and %d' %(smallerLabel, largerLabel))
                    correctionLog.write('\n\tMerging labels' + str(smallerLabel) + ' and ' + str(largerLabel))
                    if VERBOSE: print( 'Checking from label %d again.' % currentLabel )
                    correctionLog.write('\nChecking from label ' + str(currentLabel ) + ' again.')

                else:
                    if VERBOSE: print('Label ' + str(currentLabel) + ' has no large contacts.')
                    correctionLog.write('\nLabel ' + str(currentLabel) + ' has no large contacts.')
                    currentLabel = currentLabel + 1
                    if VERBOSE: print( 'Moving to next label ' + str(currentLabel) +  ' now.')
                    correctionLog.write('\nMoving to next label ' + str(currentLabel) +  ' now.')

            # Merging with area limits
            else:
                if len(largeAreaVal) != 0:
                    if VERBOSE: print('\tLabels with large area: ' + str( largeAreaLabel ) )
                    if VERBOSE: print('\tArea for above labels: ' + str( largeAreaVal ) )

                    correctionLog.write('\n\tLabels with large area: ' + str( largeAreaLabel ) )
                    correctionLog.write('\n\tArea for above labels: ' + str( largeAreaVal ) )

                    labelToMerge = largeAreaLabel[ largeAreaVal.index( min( largeAreaVal ) ) ]
                    smallerLabel = min(labelToMerge,currentLabel)
                    largerLabel = max(labelToMerge,currentLabel)
                    correctedLabelMap[np.where(correctedLabelMap == largerLabel)] = int(smallerLabel)
                    correctedLabelMap = moveLabelsUp( correctedLabelMap , largerLabel )
                    currentLabel = smallerLabel
                    lastLabel = correctedLabelMap.max()

                    if VERBOSE: print('\tMerging label %d and %d' %(smallerLabel, largerLabel))
                    correctionLog.write('\n\tMerging labels ' + str(smallerLabel) + ' and ' + str(largerLabel))
                    if VERBOSE: print( 'Checking from label %d again.' % currentLabel )
                    correctionLog.write('\nChecking from label ' + str(currentLabel ) + ' again.')

                else:
                    if VERBOSE: print('Label ' + str(currentLabel) + ' has no large contacts.')
                    correctionLog.write('\n\nLabel ' + str(currentLabel) + ' has no large contacts.')
                    currentLabel = currentLabel + 1
                    if VERBOSE: print( 'Moving to next label ' + str(currentLabel) +  ' now.')
                    correctionLog.write('\nMoving to next label ' + str(currentLabel) +  ' now.')

    if pad > 0: correctedLabelMap = removePaddingFromLabelledMap(correctedLabelMap, pad)

    timeEnd = time.time()
    timeTaken = (timeEnd - timeStart)//60

    print( '\nTime taken for correction loop: ' + str( np.round(timeTaken) ) + ' mins' )
    correctionLog.write('\nTime taken for correction loop: ' + str( np.round(timeTaken) ) + ' mins\n\n')

    correctionLog.close()

    if checkForSmallParticles == True:
        correctedCleanedLabelMap = removeSmallParticles( correctedLabelMap )
    else : correctedCleanedLabelMap = correctedLabelMap

    if saveImg == True:
        print('\nSaving corrected labelled map...')
        tiffy.imsave(outputDir+sampleName+'-correctedLabelMap.tif',correctedCleanedLabelMap.astype('uint16'))

    return correctedCleanedLabelMap

def applyPaddingToLabelledMap( labelledMap, pad ):
    paddedMap = labelledMap
    padLabMap = np.zeros( ( labelledMap.shape[0]+2*pad, labelledMap.shape[0]+2*pad, labelledMap.shape[0]+2*pad ) )
    padLabMap[pad : padLabMap.shape[0]-pad , pad : padLabMap.shape[1]-pad , pad : padLabMap.shape[ 2 ]-pad ] = labelledMap
    return padLabMap.astype(int)

def removePaddingFromLabelledMap( padLabMap, pad ):
    cleanLabMap = padLabMap[pad : padLabMap.shape[0]-pad , pad : padLabMap.shape[1]-pad , pad : padLabMap.shape[ 2 ]-pad ].astype(int)
    return cleanLabMap

def removeSmallParticles( labMapWithSmallPtcl, voxelCountThreshold = 500, saveImg=False, sampleName='', outputDir=''):
    """[TODO] can the total volume of particles be calculated - this can be used to modify the 
    gradation

    Removes particles that have a size smaller than a threshold.

    The idea is that particles that have less than 10 voxels across the
    diameter cannot be accurately measured for size and contact. Thus,
    these particles should be removed

    Assuming the diameter as the equivalent sphere diameter, the volume (voxel count)
    for such a particle will be about 500 voxels.

    Parameters
    ----------
    labMapWithSmallPtcl : ndarray
    voxelCountThreshold : unsigned integer
    saveImg : bool
    sampleName : string
    outputDir : string

    Return
    ------
    labMapUpdated : ndarray
    """
    if VERBOSE:
        print('\nRemoving small particles with voxel count smaller than ' + str( np.round( voxelCountThreshold ) ) + ' voxels' )
        print('--------------------------------------------------------------------*' )

    ptclNo = 1
    ptclCount = int( labMapWithSmallPtcl.max() )
    labMapUpdated = labMapWithSmallPtcl

    while ptclNo <= ptclCount :
        print('\nChecking ptcl no. ' + str( np.round( ptclNo ) ) )

        isolate = np.zeros_like( labMapUpdated )
        isolate[ np.where( labMapUpdated == ptclNo ) ] = 1
        voxelCount = isolate.sum()

        if voxelCount >= voxelCountThreshold :
            print('\tIts OK')
            ptclNo = ptclNo + 1

        else:
            print('\tIts smaller than threshold, removing and updating particle counts')
            labMapUpdated[ np.where( labMapUpdated  == ptclNo ) ] = int( 0 )
            labMapUpdated = moveLabelsUp( labMapUpdated , ptclNo )
            ptclCount = ptclCount - 1

    print('-------------')
    print('Complete')

    if saveImg == True:
        print('\nSaving labelled map with small particles removed')
        tiffy.imsave(outputDir+sampleName+'-noSmallCorLabMap.tif',labMapUpdated.astype('uint16'))

    return labMapUpdated

def moveLabelsUp( labelMapToFix, labelStartingWhichMoveUp ):
    """Moves the labels up - incase a label is deleted or it is merged

    This module is generally used with correction of oversegmentation or 
    removal of edge labels

    Parameters
    ----------
    labelMapToFix : ndarray
        The labe which needs to be adjusted
    labelStartingWhichMoveUp : unsigned integer
        The label such that all labels equal to or larger than this label
        are moved up in the label map

    Return
    ------
    fixedLabelMap : ndarray
        labelMap with the labels moved up after chosen label value
    """
    print( '\tUpdating Labels after %d' % labelStartingWhichMoveUp )
    deltaMatrix = np.zeros_like( labelMapToFix )
    deltaMatrix[ np.where( labelMapToFix > labelStartingWhichMoveUp ) ] = 1
    fixedLabelMap = labelMapToFix - deltaMatrix
    return fixedLabelMap

def removeEdgeLabels( labelledMapForEdgeLabelRemoval, pad=0, sampleName='', saveImg=False, outputDir='' ):
    """
    Description:
        Removes edge labels that may be cut due to the subregion boundary

    Parameters:
        labelledMapForEdgeRemoval
        pad (default set to 0)

    Return:
        labelled map with edge labels removed.
    """
    labMap = labelledMapForEdgeLabelRemoval
    numberOfLabels = labMap.max()
    startingNumberOfLabels = numberOfLabels
    currentLabel = 1

    print( '\n\nRemoving Edge Labels...' )
    print( '----------------------------*' )
    print( '\nStarting total number of labels = ' + str( labMap.max() ) + '\n' )

    startTime = time.time()
    while currentLabel <= numberOfLabels:

        if VERBOSE: print( '\n\tChecking label #' + str( currentLabel ) )
        edgeLabel = checkIfEdgeLabel( labMap,currentLabel , pad )

        if edgeLabel == True:
            if VERBOSE: print( '\tLabel #' + str( currentLabel ) + ' is an on the edge: REMOVING' )
            labMap[ np.where( labMap == currentLabel ) ] = int( 0 )
            if VERBOSE: print( '\tShifting labels up...' )
            labMap = moveLabelsUp( labMap , currentLabel )
            numberOfLabels = labMap.max()

        else:
            if VERBOSE: print( '\tLabel #' + str( currentLabel ) + ' is not on the edge: KEEPING' )
            currentLabel += 1

    print( '\nRemoval of edge labels complete' )
    print( 'Total Number of labels at start:' + str( startingNumberOfLabels ) )
    print( 'Total Number of labels remaining:' + str( numberOfLabels ) + '\n' )
    endTime = time.time()
    timeTaken = endTime - startTime
    print( 'Total time taken:~' + str( np.round( timeTaken // 60 ) ) + ' mins' )

    if saveImg == True:
        print('\nSaving corrected label map...')
        tiffy.imsave(outputDir+sampleName+'-noEdgeCorrectedLabelMap.tif', labMap.astype('uint16'))

    return labMap

def checkIfEdgeLabel( labelledMap, label, pad=0 ):
    """Checks if the particle corresponding to the label  comes in contact
    with the boundary of the scan

    The locations of the label are extracted from the labelled map. Then the x,
    y and z values of the locations are checked if they are on the lower edge
    (0) or on the upper edge of the scan

    If a padding (pad) is passed, it means that the original labelledMap has a
    padding around it and so, the lower edge is larger than zero and the upper
    edge is lower than the upper limit of the shape of the array.

    Parameters
    ----------
    labelledMap : ndarray
        Labelled map from in which the label to check exists
    label : unsigned integer
        label to check if it is located at the edge
    pad : unsigned integer
        padding that currently exists on the labelledMap

    Return
    ------
    True if the label exists on the edge and False if it does not

    """
    pointCloudArray = np.where(labelledMap == label)

    z = pointCloudArray[0].reshape( 1 , pointCloudArray[ 0 ].shape[ 0 ] )
    y = pointCloudArray[1].reshape( 1 , pointCloudArray[ 1 ].shape[ 0 ] )
    x = pointCloudArray[2].reshape( 1 , pointCloudArray[ 2 ].shape[ 0 ] )

    zLower = 0 + pad
    yLower = 0 + pad
    xLower = 0 + pad

    zUpper = labelledMap.shape[0]-1 - pad
    yUpper = labelledMap.shape[1]-1 - pad
    xUpper = labelledMap.shape[2]-1 - pad

    lowerEdge = np.any( z == zLower ) + np.any( y == yLower  ) + np.any( x == xLower  )
    upperEdge = np.any( z == zUpper ) + np.any( y == yUpper  ) + np.any( x == xUpper  )

    if lowerEdge + upperEdge > 0: return True
    else: return False

def countEdgeLabels( labelledMap ):
    countTrue = 0
    countFalse = 0
    for i in range(1,labelledMap.max() + 1):

        edge = checkIfEdgeLabel(labelledMap,i,pad=0)

        if VERBOSE: print('Label ' + str(i) + ' is edge? : ' + str(edge))

        if edge == True: countTrue+=1
        else: countFalse+=1

    print('Total Number of labels: ' + str( labelledMap.max() ) )
    print('Edge labels: ' + str( countTrue ) )
    print('Non-edge labels: ' + str( countFalse ) )

def removeLabelAndUpdate( labMap,label ):
    labMap[np.where(labMap == label)] = 0
    updatedLabMap = moveLabelsUp(labMap,label)
    return updatedLabMap
