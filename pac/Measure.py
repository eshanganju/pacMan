"""
Measure module

Written by eg
"""

import numpy as np
import math
import statistics
import scipy
import spam.label as slab
from uncertainties import unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *

from pac import Segment

# This is to plot all the text when running the functions
VERBOSE = True

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

def getParticleSize( labelledMapForParticleSizeAnalysis, calibrationFactor=1, sampleName='',saveData=False,outputDir=''):
    """
    Description:
        Computes particle size parameters for all the labels in the data

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
        #feretMax, feretMin = getMinMaxFeretDiaSPAM( labelMap=labelledMapForParticleSizeAnalysis,
        #                                            label=int(particleNum),
        #                                             numOrts=100 )
        #particleSizeDataSummary[particleNum, 6] = feretMax * calibrationFactor
        #particleSizeDataSummary[particleNum, 7] = feretMin * calibrationFactor

    if saveData == True:
        print('\nSaving particle size list...')
        np.savetxt( outputDir + sampleName + '-particleSizeList.csv',particleSizeDataSummary, delimiter=',')

    return particleSizeDataSummary  # [ Label, Volume(vx), Size0(px or mm), Size1(px or mm), Size2(px or mm), Size3(px or mm), Size4(px or mm), Size5(px or mm)]

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
    pointCloud = getZYXLocationOfLabel(labelMap,label)

    # Get covariance matrix of particle point cloud
    covarianceMatrix = np.cov( pointCloud.T )

    # Covariance matrix
    eigval, eigvec = np.linalg.eig( covarianceMatrix )

    meanZ = np.average( pointCloud[ :, 0 ] )
    meanY = np.average( pointCloud[ :, 1 ] )
    meanX = np.average( pointCloud[ :, 2 ] )

    meanMatrix = np.zeros_like( pointCloud )

    meanMatrix[ :, 0 ] = meanZ
    meanMatrix[ :, 1 ] = meanY
    meanMatrix[ :, 2 ] = meanX

    centeredLocationData = pointCloud - meanMatrix

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


def getMinMaxFeretDia( labelMap, label, numOrts=400 ):
    """
    Description:

    Parameter:

    Return:

    """
    print('Calculating Feret diameters for label ' + str(np.round(label)))
    #


def getGrainSizeDistribution( psSummary, sizeParam=1 ):
    """
    Description:

    parameter:

    Return:

    """
    print('\nGetting GSD for size param #' + str( sizeParam ) )
    label = psSummary[ : , 0 ].reshape( psSummary.shape[ 0 ] , 1 )
    vol = psSummary[ : , 1 ].reshape( psSummary.shape[ 0 ] , 1 )
    size = psSummary[: , sizeParam + 1 ].reshape( psSummary.shape[ 0 ] , 1 )
    gss = np.append( label, vol , 1 ).reshape( psSummary.shape[ 0 ] , 2 )
    gss = np.append( gss , size , 1 ).reshape( psSummary.shape[ 0 ] , 3 )

    gss = gss[ gss[ : , 2 ].argsort() ]
    totalVol = np.sum( gss[ : , 1 ] )
    pp = ( np.cumsum( gss[ : , 1 ] ) / totalVol * 100 ).reshape( gss.shape[ 0 ] , 1 )

    gsdPP = np.append( gss , pp, 1 )
    gsdPP = np.delete(gsdPP,0,0) # Removes the smallest particle (0) that comes from the getParticleSize code
    print('Done')
    return gsdPP # [ Label, Volume(vx), Size(px), percent passing(%) ]


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


def relBreak( gsdOriginal, gsdCurrent , maxSize=None , fracDim=None ):
    """
    """
    if maxSize == None: maxSize = float(input('Enter max size of particles (mm): '))
    if fracDim == None: fracDim = float(input('Enter fractal dimension to use: '))
    gsdOrig, gsdCur, gsdUlt = formatGradationsAndGetUltimate(gsdOriginal,gsdCurrent,maxSize,fracDim)
    Bp = getAreaBetweenGSDs( gsdUlt , gsdOrig )
    B = getAreaBetweenGSDs( gsdCur , gsdOrig )
    print('Bp =' + str(Bp))
    print('B = ' + str(B))
    Br = B/Bp * 100
    if VERBOSE: print('Relative breakage (Br) is: ' + str(np.round(Br,2)) + '%')
    return gsdOrig, gsdCur, gsdUlt, Br


def getAspectRatioSphericity( ):
    """
    """
    print("Measuring particle morphology...")


def getAreaBetweenGSDs( gsdUp,gsdDown,bins=1000 ):
    """
    """
    x1 = gsdUp[ : , 0 ]
    y1 = gsdUp[ : , 1 ]
    x2 = gsdDown[ : , 0 ]
    y2 = gsdDown[ : , 1 ]

    logX1 = np.log10( x1 )
    logX2 = np.log10( x2 )

    incrLogX = ( logX1.max() - logX1.min() ) / bins
    logXInterp = np.arange( logX1.min() , logX1.max() , incrLogX )
    y1Interp = np.zeros_like( logXInterp )
    y2Interp = np.zeros_like( logXInterp )
    area = np.zeros_like( logXInterp )

    for i in range( 0 , bins ):
        y1Interp[ i ] = np.interp( logXInterp[ i ] , logX1 , y1 )
        y2Interp[ i ] = np.interp( logXInterp[ i ] , logX2 , y2 )

        if i == 0:
            area[ i ] = ( y1Interp[ i ] - y2Interp[ i ] ) * incrLogX / 2

        elif i == bins - 1:
            area[ i ] = ( y1Interp[ i ] - y2Interp[ i ] ) * incrLogX / 2

        else:
            area[ i ] = ( y1Interp[ i ] - y2Interp[ i ] ) * incrLogX

    totalArea = np.sum( area )

    return totalArea


def formatGradationsAndGetUltimate( gsdOriginal,gsdCurrent,maxSize,fracDim ):
    """
    """
    xmax = float( maxSize )
    xmin = float( getStartOfUltimateGradation( xmax , fracDim ) )
    ymax = float( 100 )
    ymin = float( 0 )

    gsdCurrent = gsdCurrent.astype( float )
    gsdOriginal = gsdOriginal.astype( float )

    top = np.array( [ xmin , ymin ] ).reshape( 1 , 2 )
    bottom = np.array( [ xmax , ymax ] ).reshape( 1 , 2 )

    top = top.astype( float )
    bottom = bottom.astype( float )

    corGsdCurr = np.append( np.append( top , gsdCurrent , axis = 0 ), bottom, axis=0)
    corGsdOrig = np.append( np.append( top , gsdOriginal, axis = 0 ), bottom, axis=0)
    corGSDUlt = getUltimateGradation( xmax , xmin , fracDim )

    return corGsdOrig, corGsdCurr, corGSDUlt


def getStartOfUltimateGradation( dm , fracDim ):
    """
    """
    if VERBOSE: print( '\nGetting minimum size of ultimate gradation.' )
    d = dm
    small = False

    while small != True:
        pp = ( d / dm ) ** ( 3 - fracDim )
        if VERBOSE: print( '\tParticle size is: ' + str( np.round( d , 8 ) ) + 'mm')
        if VERBOSE: print( '\tPercentage passing is: ' + str( np.round( pp * 100 , 2 ) ) + '%\n' )
        if pp > ( 0.1 / 100 ) : d *= 0.5
        else: small = True

    return d


def getUltimateGradation( xmax , xmin , fracDim ):
    """
    """
    logIncr = float(( np.log10( xmax ) - np.log10( xmin ) ) / 100)
    logx = np.arange( np.log10( xmin ) , np.log10( xmax ) + logIncr , logIncr ).astype(float)
    x = 10**logx
    if x[-1]>xmax: x[-1]=xmax
    y = ( ( x / xmax ) ** ( 3 - fracDim ) ) * 100
    x = x.reshape( x.shape[ 0 ] , 1 )
    y = y.reshape( y.shape[ 0 ] , 1 )
    ultimateGradation = np.append( x , y , axis=1 )
    return ultimateGradation.astype( float )


def contactNormalsSpam( labelledMap, method=None ):
    """
    """
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
        print( "\tMeasuring contact using ITK watershed\n" )
        ortTabSandITK = slab.contacts.contactOrientationsAssembly( labelledData , binaryData , contactingLabels , watershed="ITK" )
        tempOrtsITK = np.zeros_like( ortTabSandITK )
        ortOnlySandITK = ortTabSandITK[ : , 2 : 5 ]

        j = 0
        for i in range( 0 , ortTabSandITK.shape[ 0 ] ):

            if ortOnlySandITK[ i , 0 ] < 0:
                ortOnlySandITK[ i ] *= -1

            if ( ortOnlySandITK[ i ] ** 2 ).sum() <= 0.999:
                print( "Contact deleted - small contact" )

            else:
                tempOrtsITK[ j ] = ortTabSandITK[ i ]
                j = j + 1

        contactTableITK = tempOrtsITK[ 0 : j , : ]

    if method == 'rw' : contTable = contactTableRW
    elif method == 'itk' : contTable = contactTableITK

    return contTable


def fabricVariablesSpam( contactTable ):
    """
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

