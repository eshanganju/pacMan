"""
Measure module

Written by eganju
"""

import numpy as np
import math
import statistics
import scipy
from scipy.ndimage import binary_erosion as erode
import spam.label as slab
from uncertainties import unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *
from pac import Segment

# This is to plot all the text when running the functions
VERBOSE = True
TESTING = False

def gsdAll( labelledMap , calib = 1.0 ):
    """
    """
    gss = getParticleSize( labelledMap ) # [ Label, Volume(vx), Size1(px), Size2(px), Size3(px), Size4(px), Size5(ND), Size6(ND)]

    gsd1 = getGrainSizeDistribution( gss , sizeParam=1 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd2 = getGrainSizeDistribution( gss , sizeParam=2 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd3 = getGrainSizeDistribution( gss , sizeParam=3 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd4 = getGrainSizeDistribution( gss , sizeParam=4 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd5 = getGrainSizeDistribution( gss , sizeParam=5 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd6 = getGrainSizeDistribution( gss , sizeParam=6 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]

    # Size in mm
    gsd1[:,2] *= calib
    gsd2[:,2] *= calib
    gsd3[:,2] *= calib
    gsd4[:,2] *= calib
    gsd5[:,2] *= calib
    gsd6[:,2] *= calib

    return gsd1[:,-2:], gsd2[:,-2:], gsd3[:,-2:], gsd4[:,-2:], gsd5[:,-2:], gsd6[:,-2:]# [ Size(mm), percent passing(%) ]


def getParticleSize( labelledMapForParticleSizeAnalysis, calibrationFactor=1,
                     sampleName='',saveData=False,outputDir=''):
    """
    Description:
        Computes particle size paameters for all the labels in the data

    Parameters:
        labelledMapForParticleSizeAnalysis,
        calibrationFactor=1,
        sampleName='',
        saveData=False,outputDir=''

    Return:
        Particle size array containing columns:
            [0] Index
            [1] Volume
            [2] Eqsp
            [3] Centroidal - max
            [4] Centroidal - med
            [5] Centroidal - min
            [6] Feret-min X
            [7] Feret-max X
    """
    numberOfParticles = int( labelledMapForParticleSizeAnalysis.max() )

    particleSizeDataSummary = np.zeros( ( numberOfParticles + 1 , 8 ) )
    print( '\nStarting measurement of particles' )
    print( '------------------------------------*' )

    for particleNum in range( 1, numberOfParticles + 1 ):
        print( "Computing size of", particleNum, "/", numberOfParticles, "particle" )
        particleSizeDataSummary[particleNum, 0] = particleNum

        # Equivalent sphere diameter
        vol, eqspDia = getEqspDia( labelMap=labelledMapForParticleSizeAnalysis, label=int(particleNum) )
        particleSizeDataSummary[particleNum, 1] = vol * (calibrationFactor**3)
        particleSizeDataSummary[particleNum, 2] = eqspDia * calibrationFactor

        # Centroidal axes lengths
        caMax, caMed, caMin = getPrincipalAxesLengths( labelMap=labelledMapForParticleSizeAnalysis,
                                                       label=int(particleNum) )
        particleSizeDataSummary[particleNum, 3] = caMax * calibrationFactor
        particleSizeDataSummary[particleNum, 4] = caMed * calibrationFactor
        particleSizeDataSummary[particleNum, 5] = caMin * calibrationFactor

        # Feret diameters
        feretMax, feretMin = getMinMaxFeretDia( labelMap=labelledMapForParticleSizeAnalysis,
                                                label=int(particleNum), numOrts=100, numRots=10)
        particleSizeDataSummary[particleNum, 6] = feretMax * calibrationFactor
        particleSizeDataSummary[particleNum, 7] = feretMin * calibrationFactor

    # Removing the extra zero in the summary (first row)
    particleSizeDataSummary = np.delete( particleSizeDataSummary, 0, 0 )

    if saveData == True:
        if VERBOSE: print('\nSaving particle size list...')
        np.savetxt( outputDir + sampleName + '-particleSizeList.csv',particleSizeDataSummary, delimiter=',')

    # [ Label, Volume(vx), Size0(px or mm), Size1(px or mm), Size2(px or mm), Size3(px or mm), Size4(px or mm), Size5(px or mm)]
    return particleSizeDataSummary


def getEqspDia( labelMap, label ):
    """
    Description:
        Fuction calculates the equivalent sphere diameter of a volume
        the equivalent sphere diameter is the diameter of the sphere with the sample volume as the partcle

    Parameter:
        labelMap
        label

    Return:
        volume of particle and the equivalent sphere diameter in pixel units
    """
    labOnlyMap = np.zeros_like( labelMap )
    labOnlyMap[np.where(labelMap == label)] = 1
    volume = labOnlyMap.sum()
    eqspLabelDia = (6*volume/(math.pi))**(1/3)

    return volume, eqspLabelDia


def getPrincipalAxesLengths( labelMap, label ):
    """
    Description:
        Computes the principal axes lengths of the particle

    Parameter:
        labelMap
        label

    Return:
       centroidalAxesLengthMax
       centroidalAxesLengthMed
       centroidalAxesLengthMin
    """
    # Centroidal axes lengths
    # Get z, y, x locations of particle
    zyxofLabel = getZYXLocationOfLabel(labelMap,label)

    # Get covariance matrix of particle point cloud
    covarianceMatrix = np.cov( zyxofLabel.T )

    # Covariance matrix
    eigval, eigvec = np.linalg.eig( covarianceMatrix )

    # centering
    meanZ, meanY, meanX = getCenterOfGravityFromZYXLocations( zyxLocationData=zyxofLabel )

    meanMatrix = np.zeros_like( zyxofLabel )

    meanMatrix[ :, 0 ] = meanZ
    meanMatrix[ :, 1 ] = meanY
    meanMatrix[ :, 2 ] = meanX

    centeredLocationData = zyxofLabel - meanMatrix

    rotationMatrix = np.zeros( ( 3, 3 ) )

    rotationMatrix[ :, 0 ] = eigvec[ 0 ]
    rotationMatrix[ :, 1 ] = eigvec[ 1 ]
    rotationMatrix[ :, 2 ] = eigvec[ 2 ]

    rotCentPointCloud = ( np.matmul( rotationMatrix, centeredLocationData.T ) ).T

    caDims = np.zeros( ( 3, 1 ) )
    caDims[ 0 ] = rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 0 ].min()
    caDims[ 1 ] = rotCentPointCloud[ :, 1 ].max() - rotCentPointCloud[ :, 1 ].min()
    caDims[ 2 ] = rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 2 ].min()

    centroidalAxesLengthMax = max( caDims )[ 0 ]
    centroidalAxesLengthMin = min( caDims )[ 0 ]
    centroidalAxesLengthMed = statistics.median( caDims )[ 0 ]
    return centroidalAxesLengthMax, centroidalAxesLengthMed, centroidalAxesLengthMin


def getMinMaxFeretDiaSPAM( labelMap, label, numOrts=100 ):
    """
    Description:

    parameter:

    Return:

    """
    labOneOnly = np.zeros_like( labelMap )
    labOneOnly[ np.where( labelMap == label ) ] = 1
    feretDims, feretOrts = slab.feretDiameters( labOneOnly, numberOfOrientations=numOrts )

    feretMax = feretDims[1,0] # Indexing start at 1 cuz spam.feret measures the 0 index - void space
    feretMin = feretDims[1,1] # ^^
    return feretMax, feretMin


def getMinMaxFeretDia( labelMap, label, numOrts=100, numRots=10 ):
    """
    Description:
        This gets the minimum and the maximim feret diameters of a particle
        The voxels of the particle are rotated along different unit vectors,
            and the length of the particles along the 3 axes are measured
        The minimum and the maximum lengths measured are returned

    Parameter:
        labelMap
        label
        numOrts: number of orientations
        numRots: number of rotation angles each orientation

    Return:
        minFeretDia
        maxFeretDia
    """
    if TESTING: print('\nChecking feret diameters of ' + str(np.round(label)) + '/' + str( np.round( labelMap.max() ) ) )

    # Get Z Y X locations of the label in an n x 3 arrangement (n is number of
    # labels in the labelMap passed)
    zyxLocations = getZYXLocationOfLabel(labelMap, label)

    # Get Z Y X locations of the center of gravity (COG) of the label
    cogZ, cogY, cogX = getCenterOfGravityFromZYXLocations ( zyxLocationData=zyxLocations )

    # Get surface voxels of the label
    surfaceOnlyLabMap = getSurfaceVoxelsofParticle(labelMap, label, returnWithLabel = True)

    # Get Z Y X locations of the surface voxels of the label in an n x 3 arrangement 
    zyxSurfaceLocations = getZYXLocationOfLabel(surfaceOnlyLabMap, label)

    # Offset the Z Y X locations of the surface voxels by the COG Z Y X
    meanMatrix = np.zeros_like(zyxSurfaceLocations)
    meanMatrix[:,0] = cogZ
    meanMatrix[:,1] = cogY
    meanMatrix[:,2] = cogX
    centeredZYXSurfaceLocations = zyxSurfaceLocations - meanMatrix

    # Make numpy array of unit vectors using Saaf and Kuiklaars method
    unitVectors = getOrientationUsingSandK( numberOfOrientations=numOrts)
    anglesInRad = np.arange( 0, math.pi+0.00001, math.pi/numRots )

    if TESTING:
        print('Orientaions in radians:')
        print(anglesInRad)

    minFeretDiameter = 0
    maxFeretDiameter = 0

    # Rotate the centered Z Y X locations of the surface voxels and measure feret diameters
    for i in range(0,unitVectors.shape[0]):
        if TESTING: print('\tChecking along orientation: ' + str(unitVectors[i]))
        a = unitVectors[ i, 0 ]
        b = unitVectors[ i, 1 ]
        c = unitVectors[ i, 2 ]
        d = ( b**2 + c**2 ) ** 0.5

        Rz = np.identity(3)
        Rz[ 1, 1 ] = c/d
        Rz[ 1, 2 ] = -b/d
        Rz[ 2, 1 ] = b/d
        Rz[ 2, 2 ] = c/d

        RzInv = np.identity(3)
        RzInv[ 1, 1 ] = c/d
        RzInv[ 1, 2 ] = b/d
        RzInv[ 2, 1 ] = -b/d
        RzInv[ 2, 2 ] = c/d

        Ry = np.identity(3)
        Ry[ 0, 0 ] = d
        Ry[ 0, 2 ] = -a
        Ry[ 2, 0 ] = a
        Ry[ 2, 2 ] = d

        RyInv = np.identity(3)
        RyInv[ 0, 0 ] = d
        RyInv[ 0, 2 ] = a
        RyInv[ 2, 0 ] = d
        RyInv[ 2, 2 ] = -a

        preRotatedZYXSurfaceLocations = np.matmul( Ry, np.matmul( Rz, centeredZYXSurfaceLocations.T ) )

        for angle in anglesInRad:
            if TESTING: print('\t\tChecking at rotation of ' + str(angle) + ' radians')
            Rx = np.identity( 3 )
            Rx[ 0, 0 ] = math.cos( angle )
            Rx[ 0, 1 ] = -math.sin( angle )
            Rx[ 1, 0 ] = math.sin( angle )
            Rx[ 1, 1 ] = math.cos( angle )

            rotatedZYXSurfaceLocations = np.matmul( Rx, preRotatedZYXSurfaceLocations )
            correctedRotatedZYXSurfaceLocations = np.matmul( RzInv,np.matmul( RyInv, rotatedZYXSurfaceLocations ) )

            sizeRange = np.zeros( ( 3, 1 ) )
            sizeRange[ 0, 0 ] = correctedRotatedZYXSurfaceLocations[ 0,: ].max() - correctedRotatedZYXSurfaceLocations[ 0, : ].min()
            sizeRange[ 1, 0 ] = correctedRotatedZYXSurfaceLocations[ 1,: ].max() - correctedRotatedZYXSurfaceLocations[ 1, : ].min()
            sizeRange[ 2, 0 ] = correctedRotatedZYXSurfaceLocations[ 2,: ].max() - correctedRotatedZYXSurfaceLocations[ 2, : ].min()

            if minFeretDiameter != 0:
                if minFeretDiameter > sizeRange.min():
                    minFeretDiameter = sizeRange.min()
                    if TESTING: print('\t\t\tUpdated minimum feret diameter to: ' + str(np.round(minFeretDiameter)))

            if maxFeretDiameter != 0:
                if maxFeretDiameter < sizeRange.max():
                    maxFeretDiameter = sizeRange.max()
                    if TESTING: print('\t\t\tUpdated max feret diameter to: ' + str(np.round(maxFeretDiameter)))

            if minFeretDiameter == 0:
                minFeretDiameter = min( sizeRange )
                if TESTING: print('\t\t\tInitial minimum feret diameter: ' + str(np.round(minFeretDiameter)))


            if maxFeretDiameter == 0:
                maxFeretDiameter = max( sizeRange )
                if TESTING: print('\t\t\tInitial max feret diameter: ' + str(np.round(maxFeretDiameter)))

    return maxFeretDiameter, minFeretDiameter


def getCenterOfGravityFromZYXLocations( zyxLocationData ):
    """

    """
    # centering
    centerOfGravityZ = np.average( zyxLocationData[ :, 0 ] )
    centerOfGravityY = np.average( zyxLocationData[ :, 1 ] )
    centerOfGravityX = np.average( zyxLocationData[ :, 2 ] )

    return centerOfGravityZ, centerOfGravityY, centerOfGravityX


def getSurfaceVoxelsofParticle( labelMap, label, returnWithLabel = True):
    """

    """
    originalParticleMap = np.zeros_like( labelMap )
    originalParticleMap[np.where(labelMap == label)] = 1
    erodedParticleMap = erode(originalParticleMap)
    surfaceMap = originalParticleMap - erodedParticleMap

    if returnWithLabel == False: return surfaceMap
    if returnWithLabel == True: return surfaceMap*int(label)


def getOrientationUsingSandK(numberOfOrientations=100):
    """
    Description:
        Picked from SPAM as is

    Parameters:

    Return:
    """
    M = int(numberOfOrientations)*2
    s = 3.6 / math.sqrt(M)
    delta_z = 2 / float(M)
    z = 1 - delta_z/2

    longitude = 0
    points = np.zeros( (numberOfOrientations,3) )

    for k in range( numberOfOrientations ):
        r = math.sqrt( 1 - z*z )
        points[k,2] = math.cos( longitude ) * r     #X
        points[k,1] = math.sin( longitude ) * r     #Y
        points[k,0] = z                             #Z
        z = z - delta_z
        longitude   = longitude + s/r

    return points


def getParticleSizeDistribution( psSummary, sizeParam='feretMin',
                              sampleName='', saveData=True, outputDir='' ):
    """Generates the particle size distribution from list of labels and sizes

    Different size parametres can be used to generate the grain size distribution.
    The size that most accurately matches the size distribution from the sieve
    analysis should ideally be used for the assessment of particle size distribution

    Parameters
    ----------
    psSummary : n by 8 np array
        This contains the results of the getParticleSize function. The array
        should have rows equal to the number of particles in the samples,
        and one column each for (0) Label, (1) Volume, (2) equivalent spere
        diameter, (3) maximum centroidal axes length, (4) intermediate
        centroidal axes length, (5) minimum centroidal axes length, (6)
        maximum feret diameter, and (7) minimum feret diameter.
    sizeParam : string
        'eqsp' - for equivalent sphere diameter
        'caMax' - for max centroidal axes length
        'caMed' - for intermediate centroidal axes length
        'caMin' - for minimum centroidal axes length
        'feretMax' - for max feret diameter
        'feretMin' - for min feret diameter
    sampleName : string
    saveData : bool
    outputDir : string

    Return
    ------
    gsdPP : n by 2 np array
        x column is the particle size and y column is the percentage passing


    """
    if VERBOSE: print( '\nGetting GSD for assuming ' + str( sizeParam ) )

    if sizeParam == 'eqsp' : sizeCol = int(2)
    elif sizeParam == 'caMax' : sizeCol = int(3)
    elif sizeParam == 'caMed' : sizeCol = int(4)
    elif sizeParam == 'caMin' : sizeCol = int(5)
    elif sizeParam == 'feretMax' : sizeCol = int(6)
    elif sizeParam == 'feretMin' : sizeCol = int(7)

    # Extracting relevant columns from particle size summary
    label = psSummary[ : , 0 ].reshape( psSummary.shape[ 0 ] , 1 )
    vol = psSummary[ : , 1 ].reshape( psSummary.shape[ 0 ] , 1 )
    size = psSummary[: , sizeCol ].reshape( psSummary.shape[ 0 ] , 1 )
    gss = np.append( label, vol , 1 ).reshape( psSummary.shape[ 0 ] , 2 )
    gss = np.append( gss , size , 1 ).reshape( psSummary.shape[ 0 ] , 3 )

    # Sorting by chosen particle size parameter
    gss = gss[ np.argsort( gss[ : , 2 ] ) ]

    # Getting percentage passing for each particle in the list
    totalVol = np.sum( gss[ : , 1 ] )
    pp = ( np.cumsum( gss[ : , 1 ] ) / totalVol * 100 ).reshape( gss.shape[ 0 ] , 1 )

    # Adding percentage passing column to the gss array
    gsdPP = np.append( gss , pp, 1 )

    if saveData == True:
        if VERBOSE: print('\nSaving particle size distribution...')
        np.savetxt( outputDir + sampleName + '-particleSizeDist.csv', gsdPP, delimiter=',')

    return gsdPP


def computeVolumeOfLabel( labelledMap, label ):
    """
    """
    labelOnlyMap = np.zeros_like(labelledMap)
    labelOnlyMap[np.where(labelledMap == label)] = 1
    volumeOfLabel = labelOnlyMap.sum()
    return volumeOfLabel


def getZYXLocationOfLabel( labelledMap, label ):
    """
    """
    particleLocationArrays = np.where(labelledMap == label)
    zyxLocationData = np.zeros( ( particleLocationArrays[ 0 ].shape[ 0 ], 3 ) )
    zyxLocationData[:,0] = particleLocationArrays[0]
    zyxLocationData[:,1] = particleLocationArrays[1]
    zyxLocationData[:,2] = particleLocationArrays[2]
    return zyxLocationData


def getRelativeBreakageHardin( psdOriginal, psdCurrent,
                               smallSizeLimit=0.075 ):
    """Computes the relative breakage parameter according to the defintion by Hardin(1985)

    This assumes that after a particle reaches a certain threshold size (sand-silt boundary
    according to Hardin (1985)), it will not break any further

    Parameters
    ----------
    psdOriginal : n by 2 numpy array
        Original particle size distribution with particle size in col 0 and percentage
        passing in col 1. the largest particle size is controlled by this gradation

    psdCurrent : n by 2 numpy array
        Current particle size distribution with particle size in col 0 and percentage
        passing in col 1. the largest particle size is controlled by this gradation

    smallSizeLimit : unsigned float
        Lower limit of integration (mm)

    Return
    ------
    potentialBreakage : unsigned float
        Area in mm2 between the ultimate and original gradation

    currentBreakage : unsigned float
        Area in mm2 between the current and original gradation

    relativeBreakage : unsigned float
        The relative breakage parameter according to Hardin (1985)
    """
    # The upper limit of integration is determined by the largest particle size
    # in the original particle size distribution
    largeSizeLimit = psdOriginal[ :,0 ].max()

    areaUnderOriginalPSD = areaUnderPSDCurve( psdOriginal,
                                              maxSize=largeSizeLimit )

    areaUnderCurrentPSD = areaUnderPSDCurve( psdCurrent,
                                             maxSize=largeSizeLimit )

    areaUnderUltimatePSD = ( math.log10( largeSizeLimit ) -
                             math.log10( smallSizeLimit ) ) * psdOriginal[ :,1 ].max()

    potentialBreakage = areaUnderUltimatePSD - areaUnderOriginalPSD
    currentBreakage = areaUnderCurrentPSD - areaUnderOriginalPSD
    relativeBreakage = currentBreakage/potentialBreakage*100

    return potentialBreakage, currentBreakage, relativeBreakage


def getRelativeBreakageEinav( gsdOriginal, gsdCurrent , fracDim=2.6,
                              smallSizeLimit=0.001 ):
    """Computes the relative breakage parameter according to the defintion of Einav (2007)

    This parameter assumes that the largest particle size does not change with crushing and
    that the crushed sand follows a fractal gradation curve

    Parameters
    ----------
    gsdOriginal : n by 2 numpy array
        Original particle size distribution with particle size in col 0 and percentage
        passing in col 1. the largest particle size is controlled by this gradation

    gsdCurrent : n by 2 numpy array
        Current particle size distribution with particle size in col 0 and percentage
        passing in col 1.

    fracDim : unsigned float
        This is the fractal dimension used to compute the ultimate gradation

    smallSizeLimit : unsigned float
        Lower limit of integration (mm)

    Return
    ------
    potentialBreakage : unsigned float
        Area in mm2 between the ultimate and original gradation

    currentBreakage : unsigned float
        Area in mm2 between the current and original gradation

    relativeBreakage : unsigned float
        The relative breakage parameter according to Einav (2007)
    """
    # The upper limit of integration is determined by the largest particle size
    # in the original particle size distribution
    largeSizeLimit = psdOriginal[ :,0 ].max()

    areaUnderOriginalPSD = areaUnderPSDCurve( psdOriginal,
                                              maxSize=largeSizeLimit )

    areaUnderCurrentPSD = areaUnderPSDCurve( psdCurrent,
                                             maxSize=largeSizeLimit )

    psdUltimate = getUltimateFractalParticleSizeDistribution( minSize=smallSizeLimit,
                                                              maxSize=largeSizeLimit,
                                                              fractalDimension=fracDim )

    areaUnderUltimatePSD = areaUnderPSDCurve( psdUltimate,
                                              maxSize=largeSizeLimit)

    potentialBreakage = areaUnderUltimatePSD - areaUnderOriginalPSD
    currentBreakage = areaUnderCurrentPSD - areaUnderOriginalPSD
    relativeBreakage = currentBreakage/potentialBreakage*100

    return potentialBreakage, currentBreakage, relativeBreakage

def getUltimateFractalParticleSizeDistribution(minSize=0.0,maxSize=0.0,
                                               fractalDimension=2.6,num=100):
    """Gets the ultimate particle size distribution

    The ultimate particle size distribution is obtained following the fractal
    distribution

    The equation is:
        P = (d/dMax)^(3-Xi)
            where p is the percentage by mass smaller than particle size d
            dMax is the largest particle size in the original sample
            Xi is the fractal dimension of the PSD in the ultimate state of crushing

    Parameters
    ----------
    minSize : unsigned float
        smallest particle size in the PSD

    maxSize : unsigned float
        largest particle size in the PSD

    fractalDimension : unsigned float
        the fractal dimension of the ultimate gradation

    num : integer
        Number of points in the PSD

    Return
    ------
    ultimatePSDFractal : unsigned float array
        Ultimate particle size distribution with particle size in col 0 and percentage
        passing in col 1. Col 1 is ascending order.
    """
    if minSize == 0.0:
        minSize = float( input( 'Enter the min particle size (mm):' ) )
    if maxSize == 0.0:
        maxSize = float( input( 'Enter the max particle size (mm):' ) )

    ultimatePSDFractal = np.zeros((num,2))
    ultimatePSDFractal[ :, 0 ] = np.linspace( minSize, maxSize, num )
    ultimatePSDFractal[ :, 1 ] = ( ultimatePSDFractal[ :, 0 ] / maxSize ) ** ( 3-fractalDimension )

    return ultimatePSDFractal


def areaUnderPSDCurve( psd, maxSize=0.0 )
    """Computes the area under the particle size distribution curve passed to it

    Since the GSD curves are plotted on a semi-log graph, the log10 of the x
    axis (particle size) is taken before computing the area.

    Parameters
    ----------
    psd : array of floats
        Particle size distribution containing particle size (mm) in col 1 and
        percentage passing (%) in col 2. Col 1 should be in ascending order

    maxSize : unsigned float
        upper limit of integration in mm

    Return
    ------
    areaUnderPSDCurve : unsigned float
        area under the psd
    """
    areaUnderCurve = 0.0
    numberOfPoints = psd.shape[0]

    if maxSize == 0.0: maxSize = float(input('Enter the max particle size (mm): '))

    # Adding the max size to the array, min is not needed
    psd = np.append( psd, np.array( [maxSize , 100.0] ).reshape( 1, 2 ), 0 )

    for i in range( 0, numberOfPoints - 1 ):
        Y = psd[ i, 1 ]
        deltaX = np.log10( psd[ i+1, 0 ] ) - np.log10( psd[ i, 0 ] )
        deltaY = psd[ i+1, 1 ] - psd[ i, 1 ]
        areaRectangle = Y * deltaX
        areaTriangle = 0.5 * ( deltaX * deltaY )
        areaTotal = areaRectangle + areaTriangle
        areaUnderCurve = areaUnderCurve + areaTotal

    return areaUnderPSDCurve


def contactNormalsSPAM( labelledMap, method=None ):
    """This computes the orientations of the inter particle contacts
    between particles in the scanned data

    Parameters
    ----------

    Return
    ------

    """
    if VERBOSE:
        print( "\nMeasuring contact normals using SPAM library\n" )

    labelledData = labelledMap
    binaryData = np.zeros_like( labelledData )
    binaryData[ np.where( labelledData != 0 ) ] = int( 1 )

    contactVolume, Z, contactsTable, contactingLabels = slab.labelledContacts( labelledData )

    if method == None:
        method = input('Enter the method to use (itk or rw): ')

    if method == 'rw':
        print("\tMeasuring contact using Random Walker\n")
        ortTabSandRW = slab.contacts.contactOrientationsAssembly( labelledData , binaryData , contactingLabels , watershed = "RW" )
        tempOrtsRW = np.zeros_like( ortTabSandRW )
        ortOnlySandRW = ortTabSandRW[ : , 2 : 5 ]

        j = 0
        for i in range( 0 , ortTabSandRW.shape[ 0 ] ):

            if ortOnlySandRW[ i , 0 ] < 0:
                ortOnlySandRW[ i ] *= -1

            if ( ortOnlySandRW[ i ] ** 2 ).sum() <= 0.999:
                if VERBOSE: print( "Contact deleted - small contact" )

            else:
                tempOrtsRW[ j ] = ortTabSandRW[ i ]
                j = j + 1

        contactTableRW = tempOrtsRW[ 0 : j , : ]

    else :
        if VERBOSE:
            print( "\tMeasuring contact using ITK watershed\n" )
        ortTabSandITK = slab.contacts.contactOrientationsAssembly( labelledData , binaryData , contactingLabels , watershed="ITK" )
        tempOrtsITK = np.zeros_like( ortTabSandITK )
        ortOnlySandITK = ortTabSandITK[ : , 2 : 5 ]

        j = 0
        for i in range( 0 , ortTabSandITK.shape[ 0 ] ):

            if ortOnlySandITK[ i , 0 ] < 0:
                ortOnlySandITK[ i ] *= -1

            if ( ortOnlySandITK[ i ] ** 2 ).sum() <= 0.999:
                if VERBOSE:
                    print( "Contact deleted - small contact" )

            else:
                tempOrtsITK[ j ] = ortTabSandITK[ i ]
                j = j + 1

        contactTableITK = tempOrtsITK[ 0 : j , : ]

    if method == 'rw' : contTable = contactTableRW
    elif method == 'itk' : contTable = contactTableITK

    return contTable


def fabricVariablesSpam( contactTable ):
    """Computes the fabric tensors using the SPAM library

    Parameter
    ---------

    Return
    ------

    """
    orts = contactTable[ :, 2:5]
    F1, F2, F3 = slab.fabricTensor( orts )
    return F1, F2, F3


def fabricVariablesWithUncertainity( contactTable, vectUncert = 0 ):
    """
    """
    vectors = contactTable[ :, 2:5]                         # contactTable col 0 is first particle label and col 1 is second particle number
    uncertVectors = vectUncert*(np.ones_like(vectors))      # uncertainity in the vector components

    uncertVectorArray = unp.uarray(vectors,uncertVectors)   # the uncertainity libraries takes vectors and the undertainity of those vectors to make an array

    N = np.zeros((3,3))
    F = np.zeros((3,3))
    Fq = np.zeros((1,1))

    uN = unp.uarray(N,N)
    uF = unp.uarray(F,F)
    uFq = unp.uarray(Fq,Fq)

    for i in range(0,uncertVectorArray.shape[0]):
        uN[0,0] = uN[0,0] + (uncertVectorArray[i,0])*(uncertVectorArray[i,0])
        uN[0,1] = uN[0,1] + (uncertVectorArray[i,0])*(uncertVectorArray[i,1])
        uN[0,2] = uN[0,2] + (uncertVectorArray[i,0])*(uncertVectorArray[i,2])
        uN[1,0] = uN[1,0] + (uncertVectorArray[i,1])*(uncertVectorArray[i,0])
        uN[1,1] = uN[1,1] + (uncertVectorArray[i,1])*(uncertVectorArray[i,1])
        uN[1,2] = uN[1,2] + (uncertVectorArray[i,1])*(uncertVectorArray[i,2])
        uN[2,0] = uN[2,0] + (uncertVectorArray[i,2])*(uncertVectorArray[i,0])
        uN[2,1] = uN[2,1] + (uncertVectorArray[i,2])*(uncertVectorArray[i,1])
        uN[2,2] = uN[2,2] + (uncertVectorArray[i,2])*(uncertVectorArray[i,2])

    uN = uN / uncertVectorArray.shape[0]

    utraceF = uN[0,0] + uN[1,1] + uN[2,2]

    uF[0,0] = 15/2 * ( uN[0,0] - (1/3)*utraceF )
    uF[0,1] = 15/2 * ( uN[0,1] )
    uF[0,2] = 15/2 * ( uN[0,2] )
    uF[1,0] = 15/2 * ( uN[1,0] )
    uF[1,1] = 15/2 * ( uN[1,1] - (1/3)*utraceF )
    uF[1,2] = 15/2 * ( uN[1,2] )
    uF[2,0] = 15/2 * ( uN[2,0] )
    uF[2,1] = 15/2 * ( uN[2,1] )
    uF[2,2] = 15/2 * ( uN[2,2] - (1/3)*utraceF )

    uFq[0,0] = ((3/2)*( uF[0,0]*uF[0,0] + uF[0,1]*uF[0,1] + uF[0,2]*uF[0,2] + uF[1,0]*uF[1,0] + uF[1,1]*uF[1,1] + uF[1,2]*uF[1,2] + uF[2,0]*uF[2,0] + uF[2,1]*uF[2,1] + uF[2,2]*uF[2,2])) ** 0.5

    return uN, uF, uFq


def getCoordinationNumberList(labelledMap):
    """
    """
    numberOfLabels = labelledMap.max()
    coordinationNumberArray = np.zeros( ( numberOfLabels, 2) )

    labelledMap = Segment.applyPaddingToLabelledMap(labelledMap, 2)

    for currentLabel in range(2, numberOfLabels + 1):
        print('\nChecking for label ' + str(np.round(currentLabel)))
        contactLabels = slab.contactingLabels( labelledMap, currentLabel, areas=False)
        numberOfContacts = len(contactLabels)
        coordinationNumberArray[currentLabel-1,0] = currentLabel
        coordinationNumberArray[currentLabel-1,1] = numberOfContacts

    # labelledMap = Segment.removePaddingFromLabelledMap(padLabMap, 2)

    return coordinationNumberArray


