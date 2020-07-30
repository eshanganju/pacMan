'''
Measure module

'''

import numpy as np
import math
import statistics
import scipy
import spam.label as slab
from uncertainties import unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *

VERBOSE = True

def gsd( labelledMap , calib = 0.01193 ):
    gss = getParticleSize( labelledMap ) # [ Label, Volume(vx), Size1(px), Size2(px), Size3(px), Size4(px), Size5(ND), Size6(ND)]

    gsd1 = getGrainSizeDistribution( gss , sizeParam=1 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd2 = getGrainSizeDistribution( gss , sizeParam=2 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd3 = getGrainSizeDistribution( gss , sizeParam=3 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]
    gsd4 = getGrainSizeDistribution( gss , sizeParam=4 ) # [ Lable, Volume(vx), Size(px), percent passing(%) ]

    # Size in mm
    gsd1[:,2] *= calib
    gsd2[:,2] *= calib
    gsd3[:,2] *= calib
    gsd4[:,2] *= calib

    # Volume in mm3
    #gsd1[:,1] *= ( calib ** 3 )
    #gsd2[:,1] *= ( calib ** 3 )
    #gsd3[:,1] *= ( calib ** 3 )
    #gsd4[:,1] *= ( calib ** 3 )

    return gsd1[:,-2:], gsd2[:,-2:], gsd3[:,-2:], gsd4[:,-2:] # [ Size(mm), percent passing(%) ]

# TODO:Split getParticleSize in to smaller units
#   One for equivlent sphere
#   Three for centroidal axis dimensions

def getParticleSize( labelledMapForParticleSizeAnalysis ):
    numberOfParticles = int( labelledMapForParticleSizeAnalysis.max() )

    # Particle size summary columns
    # [0] Index, [1] Volume, [2] Eqsp, [3] Centroidal - max, [4] Centroidal - med, [5] Centroidal - min, [6] and [7] are open
    particleSizeDataSummary = np.zeros( ( numberOfParticles + 1 , 8 ) )
    print( "Starting measurement of particles..." )

    for particleNum in range( 1, numberOfParticles + 1 ):
        print( "Computing size of", particleNum, "/", numberOfParticles, "particle" )
        particleSizeDataSummary[particleNum, 0] = particleNum

        # Equivalent sphere diameter (Param 1)
        vol = computeVolumeOfLabel(labelledMapForParticleSizeAnalysis,particleNum)
        particleSizeDataSummary[particleNum, 1] = vol
        sphDia = 2 * ( ( ( 3 * vol ) / ( 4 * math.pi ) ) ** ( 1 / 3 ) )
        particleSizeDataSummary[particleNum, 2] = sphDia

        # Feret diameters
        # Get z, y, x locations of particle
        pointCloud = getZYXLocationOfLabel(labelledMapForParticleSizeAnalysis,particleNum)

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

        caMax = max( caDims )[ 0 ] # Param 2
        caMin = min( caDims )[ 0 ] # Param 3
        caMed = statistics.median( caDims )[ 0 ] # Param 4

        particleSizeDataSummary[particleNum, 3] = caMax
        particleSizeDataSummary[particleNum, 4] = caMed
        particleSizeDataSummary[particleNum, 5] = caMin

    return particleSizeDataSummary  # [ Label, Volume(vx), Size1(px), Size2(px), Size3(px), Size4(px), Size5(ND), Size6(ND)]

def getGrainSizeDistribution(psSummary,sizeParam=1):
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

def computeVolumeOfLabel( labelledMap, label):
    labelOnlyMap = np.zeros_like(labelledMap)
    labelOnlyMap[np.where(labelledMap == label)] = 1
    volumeOfLabel = labelOnlyMap.sum()
    return volumeOfLabel

def getZYXLocationOfLabel( labelledMap, label):
    particleLocationArrays = np.where(labelledMap == label)
    zyxLocationData = np.zeros( ( particleLocationArrays[ 0 ].shape[ 0 ], 3 ) )
    zyxLocationData[:,0] = particleLocationArrays[0]
    zyxLocationData[:,1] = particleLocationArrays[1]
    zyxLocationData[:,2] = particleLocationArrays[2]
    return zyxLocationData

def relBreak( gsdOriginal, gsdCurrent , maxSize=None , fracDim=None ):
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

def getAspectRatioSphericity( particleSizeSummary ):
    print("Measuring particle morphology...")
    '''
    Take the ratio of the long axis - could be the ratio of CA max and CA min
    '''

def getAreaBetweenGSDs(gsdUp,gsdDown,bins=1000):
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

def formatGradationsAndGetUltimate(gsdOriginal,gsdCurrent,maxSize,fracDim):
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
    logIncr = float(( np.log10( xmax ) - np.log10( xmin ) ) / 100)
    logx = np.arange( np.log10( xmin ) , np.log10( xmax ) + logIncr , logIncr ).astype(float)
    x = 10**logx
    if x[-1]>xmax: x[-1]=xmax
    y = ( ( x / xmax ) ** ( 3 - fracDim ) ) * 100
    x = x.reshape( x.shape[ 0 ] , 1 )
    y = y.reshape( y.shape[ 0 ] , 1 )
    ultimateGradation = np.append( x , y , axis=1 )
    return ultimateGradation.astype( float )

def contactNormalsSpam( labelledMap, method=None):
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

def fabricVariables( contactTable ):
    orts = contactTable[ :, 2:5]
    F1, F2, F3 = slab.fabricTensor( orts )
    return F1, F2, F3

def fabricVariablesWithUncertainity( contactTable, vectUncert = 0.26 ):
    vectors = contactTable[ :, 2:5]
    uncertVectors = vectUncert*(np.ones_like(vectors))

    uncertVectorArray = unp.uarray(vectors,uncertVectors)

    N = np.zeros((3,3))
    F = np.zeros((3,3))
    Fq = 0.0

    uN = unp.uarray(N,N)
    uF = unp.uarray(F,F)
    uFq = ufloat(Fq,Fq)

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

    utraceF = N[0,0] + N[1,1] + N[2,2]


    uF[0,0] = 15/2 * ( uN[0,0] - (1/3)*utraceF )
    uF[0,1] = 15/2 * ( uN[0,1] )
    uF[0,2] = 15/2 * ( uN[0,2] )
    uF[1,0] = 15/2 * ( uN[1,0] )
    uF[1,1] = 15/2 * ( uN[1,1] - (1/3)*utraceF )
    uF[1,2] = 15/2 * ( uN[1,2] )
    uF[2,0] = 15/2 * ( uN[2,0] )
    uF[2,1] = 15/2 * ( uN[2,1] )
    uF[2,2] = 15/2 * ( uN[2,2] - (1/3)*utraceF )

    uFq = ((3/2)*( F[0,0]*F[0,0] + F[0,1]*F[0,1] + F[0,2]*F[0,2] + F[1,0]*F[1,0] + F[1,1]*F[1,1] + F[1,2]*F[1,2] + F[2,0]*F[2,0] + F[2,1]*F[2,1] + F[2,2]*F[2,2])) ** 0.5

    return uN, uF, uFq
