'''
Segment
'''

# Importing libraries
from scipy.ndimage.morphology import distance_transform_edt as edt
from scipy.ndimage.morphology import binary_fill_holes
from scipy.ndimage.morphology import binary_opening
from skimage.morphology import local_maxima as localMaxima
from skimage.morphology import h_maxima as hmax
from skimage.morphology import watershed as wsd
from skimage.filters import threshold_otsu
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt
import numpy as np
import spam.label as slab
import time

# This is to plot all the text in the methods
VERBOSE = True


def performEDTWS( filteredGLIMap, currentVoidRatio, outputFilesLocation, sampleName):
    print('Starting EDT-WS')
    print('----------------------*\n')

    # Binarize
    print('Which binarization method to follow?')
    binarizationMethodToFollow = input('(1) Otsu threshold, (2) User input based threshold, ([3]) Density based threshold: ')

    # Otsu threshold
    if binarizationMethodToFollow == '1':
        binMap, gliThreshold = binarizeAccordingToOtsu( filteredGLIMap, outputFilesLocation, sampleName)
        voidRatioCT = calcVoidRatio( binMap )

    # User input based threshold
    elif binarizationMethodToFollow == '2':
        binMap, gliThreshold = binarizeAccordingToUserThreshold( filteredGLIMap, outputFilesLocation, sampleName )
        voidRatioCT = calcVoidRatio( binMap )

    # Density based threshold
    else:
        binMap, gliThreshold = binarizeAccordingToDensity(filteredGLIMap, currentVoidRatio, outputFilesLocation, sampleName )
        voidRatioCT = calcVoidRatio( binMap )

    # EDM
    edMap = obtainEuclidDistanceMap( binMap )

    # Markers
    edPeakMrkrMap = obtainLocalMaximaMarkers( edMap )

    # Labelled Map
    labelledMap = obtainLabelledMapUsingWaterShedAlgorithm( binMap, edMap, edPeakMrkrMap, outputFilesLocation )

    # Correction of labelled map
    lblCorrectionMethod, correctedLabelledMap = fixErrorsInSegmentation( labelledMap )

    # Returns
    return gliThreshold, binMap, voidRatioCT, edMap, edPeakMrkrMap, lblCorrectionMethod, correctedLabelledMap

def binarizeAccordingToOtsu( gliMapToBinarize ):
    print('\nRunning Otsu Binarization')
    print('----------------------------*')
    otsuThreshold = threshold_otsu( gliMapToBinarize )
    binaryMap = np.zeros_like( gliMapToBinarize )

    binaryMap[ np.where( gliMapToBinarize > otsuThreshold ) ] = 1
    binaryMap = binaryMap.astype( int )
    e1 = calcVoidRatio( binaryMap )

    binaryMap2 = fillHoles( binaryMap )
    e2 = calcVoidRatio( binaryMap2 )

    binaryMap3 = removeSpecks( binaryMap )
    e3 = calcVoidRatio( binaryMap3 )

    print( 'Global Otsu threshold = %f' % otsuThreshold )
    print( 'Void ratio after threshold = %f' % e1 )
    print( 'Void ratio after filling holes = %f' % e2 )
    print( 'Void ratio after removing specks = %f' % e3 )

    return binaryMap3

def binarizeAccordingToUserThreshold( gliMapToBinarize, outputLocation, sampleName):
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

    print( 'Global User threshold = %f' % userThreshold )
    print( 'Void ratio after threshold = %f' % e1 )  
    print( 'Void ratio after filling holes = %f' % e2 )   
    print( 'Void ratio after removing specks = %f' % e3 )    

    # Saving files
    userThresholdFileName = outputLocation + sampleName + '-userThresholdDetails.txt'
    f = open( userThresholdFileName, 'w+' )        
    f.write( 'Global user threshold = %f' % userThreshold )
    f.write( '\nVoid ratio after threshold = %f' % e1 )  
    f.write( '\nVoid ratio after filling holes = %f' % e2 )   
    f.write( '\nVoid ratio after removing specks = %f' % e3 )
    f.close() 

    print('Output file saved as ' + userThresholdFileName)
    print('---------------------*')

    return binaryMap, userThreshold

def binarizeAccordingToDensity(gliMapToBinarize, knownVoidRatio, outputLocation, sampleName):
    otsuBinMap, otsuThreshold = binarizeAccordingToOtsu(gliMapToBinarize, outputLocation, sampleName)
    currentThreshold = otsuThreshold
    currentVoidRatio = calcVoidRatio( otsuBinMap )
    targetVoidRatio = knownVoidRatio
    greyLvlMap = gliMapToBinarize

    tolerance = 0.0001
    deltaVoidRatio = targetVoidRatio - currentVoidRatio

    absDeltaThreshold = 0
    maxAbsDeltaThreshold = 500

    iterationNum = 1
    maxIterations = 50

    userThresholdStepsTextFileName = outputLocation +  sampleName + '-userThresholdDensitySteps.txt'
    f = open( userThresholdStepsTextFileName, "w+")
    f.write( "Steps of density based user threshold\n\n" )
    f.close()

    while( abs( deltaVoidRatio ) > tolerance and iterationNum <= maxIterations ):
        incrementSign = deltaVoidRatio/abs(deltaVoidRatio)  

        if abs( deltaVoidRatio ) > 0.100:
            absDeltaThreshold = 500

        elif abs( deltaVoidRatio ) > 0.05:
            absDeltaThreshold = 250

        elif abs( deltaVoidRatio ) > 0.01:
            absDeltaThreshold = 100

        elif abs( deltaVoidRatio ) > 0.005:
            absDeltaThreshold = 50

        elif abs( deltaVoidRatio ) > 0.001:
            absDeltaThreshold = 25

        elif abs( deltaVoidRatio ) > 0.0005:
            absDeltaThreshold = 10

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

        f = open( userThresholdStepsTextFileName, 'a' )
        f.write( "\nIteration %d:\n" % iterationNum)
        f.write( "Current threshold: %d\n" % round( currentThreshold ) )
        f.write( "Otsu Threshold: %d\n" %  round( otsuThreshold ) )
        f.write( "Current void ratio:  %0.5f\n" % round( currentVoidRatio,5 ) )
        f.write( "Target void ratio: %0.5f\n" % round( targetVoidRatio,5 ) )
        f.write( "Target - Current void ratio: %0.5f\n" % round( deltaVoidRatio, 5 ) )
        f.close()

        iterationNum = iterationNum + 1

    densityBasedThresholdTextFileName = outputLocation +  sampleName + '-userThresholdDensityBased.txt'
    f = open( densityBasedThresholdTextFileName, "w+" )               
    f.write( "Global density based threshold = %f\n" % currentThreshold )
    f.close()

    return currentBinaryMap

def calcVoidRatio(binaryMapforVoidRatioCalc):
    '''
    Parameters
    ----------
    binaryMap : numpy array with 1 and 0

    Returns
    -------
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
    '''
    Parameters
    ----------
    oldBinaryMap : Numpy array of binary data
        obtained from binarization

    Returns
    -------
    newBinaryMap holes in the map is filled

    '''
    newNoHoleBinaryMap = binary_fill_holes( oldBinaryMapWithHoles )
    return newNoHoleBinaryMap.astype(int) 

def removeSpecks( oldBinaryMapWithSpecks ):
    '''
    Parameters
    ----------
    oldBinaryMap : Numpy array of binary data
        obtained from binarization

    Returns
    -------
    newBinaryMap with white specks in the map removed


    '''
    newNoSpekBinaryMap = binary_opening( oldBinaryMapWithSpecks )   
    return newNoSpekBinaryMap.astype(int)

def obtainEuclidDistanceMap(binaryMapForEDM ):
    print('\nFinding Euclidian distance map (EDM)')
    print('------------------------------------*')
    edMap = edt( binaryMapForEDM )             
    print( "EDM Created" )        
    return edMap

def obtainLocalMaximaMarkers( edMapForPeaks ):
    '''
    Fix this to be better at obtaining markers
    Expected number of particles?
    Particle size?
    '''
    print('\nObtaining peaks of EDM...')
    print('------------------------------*')
    h = int( input( 'Enter the minimum height for a peak (px): ') )
    print( 'Finding local maxima in EDM' )

    edmPeakMarkers = hmax( edMapForPeaks, h ).astype(int)

    print( 'Found local maximas in EDM (greater than %i px)' % h )

    print( 'Resetting count of peaks' )
    count=0
    for frame in range( 0, edmPeakMarkers.shape[ 0 ] ):
        for row in range( 0, edmPeakMarkers.shape[ 1 ] ):
            for col in range( 0, edmPeakMarkers.shape[ 2 ] ):
                if edmPeakMarkers[ frame ][ row ][ col ] == 1:
                    edmPeakMarkers[ frame ][ row ][ col ] = count + 1
                    count = count + 1
        print( 'Processed ' + str(frame + 1) + ' out of ' + str(edmPeakMarkers.shape[ 0 ]) + ' slices' )

    return edmPeakMarkers

def obtainLabelledMapUsingITKWS( gliMap ):
    print( '\nSegmenting particles by ITK topological watershed' )
    binMask = binarizeAccordingToOtsu( gliMap )
    edMap = obtainEuclidDistanceMap( binMask )
    edPeaksMap = obtainLocalMaximaMarkers(edMap)
    labelledMap = wsd( -edMap, markers = edPeaksMap, mask = binMask )
    print( 'Watershed segmentation complete' )
    return labelledMap

def fixErrorsInSegmentation(labelledMapForOSCorr, pad=2):
    print('\n---------------------------*')
    print('Entering label correction')
    print('---------------------------*')

    areaLimit = int(input('Input area limit (px): '))

    if pad > 0: labelledMapForOSCorr = applyPaddingToLabelledMap(labelledMapForOSCorr, pad)
    lastLabel = labelledMapForOSCorr.max()
    currentLabel = 1
    correctedLabelMap = labelledMapForOSCorr
    timeStart = time.time()

    print('Currently padding is ' + str(pad) + ' px')

    considerEdgeLabels = input("Consolidate edge labels also ([y]/n): ")
    if considerEdgeLabels == 'n':
        print('\tOk, not consolidating edge labels')
        edgePad = pad
    else:
        print('\tOk, consolidating edge labels')
        edgePad = 0

    # Loop through labels till all OS labels are merged
    while currentLabel <= lastLabel:
        isEdgeLabel = checkIfEdgeLabel( correctedLabelMap, currentLabel, edgePad )

        # If edge label, move to next label
        if isEdgeLabel == True:
            print('Label %d is an edge label, moving to next label' % currentLabel)
            currentLabel += 1

        # If not edge label, check contacting labels
        else:
            contactLabel, contactArea = slab.contactingLabels(correctedLabelMap, currentLabel, areas=True)
            print('\nLabel ' + str( currentLabel ) + ' is contacting ' + str( contactLabel ) )
            print('\tWith areas: ' + str( contactArea ) )

            '''
            TODO:
                Read only those labels that have area larger than the limit - speed up the process
                When consolidating the edge lables, account for non-contact labels
            '''

            # Check if any contact is greater than threshold area
            '''
            Check if any contact area is larger than a limit
            find the first contact that is larger and merge
            '''
            for positionNumber in range(0, contactArea.shape[0]):
                if contactArea[positionNumber] > areaLimit:
                    print('Area between label %d and %d is greater than limit' %(currentLabel, contactLabel[positionNumber]))

                    if currentLabel < contactLabel[positionNumber]:
                        correctedLabelMap[np.where(correctedLabelMap == contactLabel[positionNumber])] = int(currentLabel)
                        if VERBOSE: print('\tMerging label %d and %d' %(currentLabel, contactLabel[positionNumber]))
                        correctedLabelMap = moveLabelsUp(correctedLabelMap,contactLabel[positionNumber])
                        lastLabel = correctedLabelMap.max()
                        print( 'Checking from label %d again' % currentLabel )
                        break

                    else:
                        correctedLabelMap[np.where(correctedLabelMap == currentLabel)] = int(contactLabel[positionNumber])
                        if VERBOSE: print('\tMerging label %d and %d' %(currentLabel, contactLabel[positionNumber]))
                        correctedLabelMap = moveLabelsUp(correctedLabelMap,currentLabel)
                        lastLabel = correctedLabelMap.max()
                        currentLabel = contactLabel[positionNumber]
                        print( 'Checking from label %d onwards now' % currentLabel )
                        break
                else:
                    print('Area between label %d and %d is ok' % ( currentLabel, contactLabel[positionNumber]))
                    if positionNumber == contactArea.shape[0]-1:
                        currentLabel = currentLabel + 1
                        print( 'Moving to next label %d now' % currentLabel )

    print('Label edition completed.')

    if pad > 0: correctedLabelMap = removePaddingFromLabelledMap(correctedLabelMap, pad)

    timeEnd = time.time()
    timeTaken = (timeEnd - timeStart)//60

    print( 'Time taken for correction loop: ' + str( timeTaken ) + ' mins' )

    return correctedLabelMap

def applyPaddingToLabelledMap(labelledMap, pad):
    '''
    Some spam function have problems with edge labels
    This function applys a padding of zeros (2px) to the edges 
    '''
    paddedMap = labelledMap
    padLabMap = np.zeros( ( labelledMap.shape[0]+2*pad, labelledMap.shape[0]+2*pad, labelledMap.shape[0]+2*pad ) )
    padLabMap[pad : padLabMap.shape[0]-pad , pad : padLabMap.shape[1]-pad , pad : padLabMap.shape[ 2 ]-pad ] = labelledMap
    return padLabMap.astype(int)

def removePaddingFromLabelledMap(padLabMap, pad):
    '''
    This removes padding around the 
    '''
    cleanLabMap = padLabMap[pad : padLabMap.shape[0]-pad , pad : padLabMap.shape[1]-pad , pad : padLabMap.shape[ 2 ]-pad ].astype(int)
    return cleanLabMap

def getTableOfContactAndArea(labelledMapForContactAndArea):
    '''
    RETURN:
        np array of shape: Number of contacts X 7
            contactingLabel1
            contactingLabel2
            contactArea
            normalZZ
            normalYY
            normalXX
            note
    '''

def getTableOfContactingLabels(labelledMapForContactDetection):
    '''
    RETURN:
        np Array of shape: Number of contacts X 2
            contactingLabel1
            contactingLabel2
    '''

def getContactAreaAndNormals(labelledMap, contactList):
    '''
    RETURN:
        np Array of shape Number of contact X 3
            contactingLabel1
            contactingLabel2
            contactArea
            contactNormalZZ
            contactNormalYY
            contactNormalXX
    '''

def moveLabelsUp(labelMapToFix, labelStartingWhichMoveUp):
    print('\tUpdating Labels after %d' % labelStartingWhichMoveUp)
    deltaMatrix = np.zeros_like(labelMapToFix)
    deltaMatrix[np.where(labelMapToFix > labelStartingWhichMoveUp)]=1
    fixedLabelMap = labelMapToFix - deltaMatrix
    return fixedLabelMap

def removeEdgeLabels(labelledMapForEdgeLabelRemoval,pad=0):
    labMap = labelledMapForEdgeLabelRemoval
    numberOfLabels = labMap.max()
    currentLabel = 1
    print('\n\nRemoving Edge Labels...')
    print('Starting total number of labels = ' + str(labMap.max()) + '\n')

    while currentLabel <= numberOfLabels:
        print('\nChecking label #' + str(currentLabel))
        edgeLabel = checkIfEdgeLabel(labMap,currentLabel,pad)
        if edgeLabel == True:
            if VERBOSE: print('\tLabel #' + str(currentLabel) + ' is an on the edge: REMOVING')
            labMap[np.where(labMap == currentLabel)] = int(0)
            if VERBOSE: print('\tShifting labels up...')
            labMap = moveLabelsUp(labMap,currentLabel)
            numberOfLabels = labMap.max()
        else:
            if VERBOSE: print('\tLabel #' + str(currentLabel) + ' is not on the edge: KEEP')
            currentLabel += 1

    print( 'Total Number of labels remaining:' + str( numberOfLabels ) + '\n' )
    return labMap

def checkIfEdgeLabel(labelledMap, label, pad):
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

    if lowerEdge + upperEdge > 0:
        return True
    else: return False

def countEdgeLabels(labelledMap):
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


