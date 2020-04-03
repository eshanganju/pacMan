'''
Measure class

'''

import numpy as np
import math
import statistics
import scipy
import spam.label as slab

def gsd( labelledMap , calib = 0.01193 ):
    gss = getParticleSize( labelledMap ) # [ Label, Volume(px3), Size1(px), Size2(px), Size3(px), Size4(px), Size5(ND), Size6(ND)]

    gsd1 = getGrainSizeDistribution( gss , sizeParam=1 ) # [ Lable, Volume(px3), Size(px), percent passing(%) ]
    gsd2 = getGrainSizeDistribution( gss , sizeParam=2 ) # [ Lable, Volume(px3), Size(px), percent passing(%) ]
    gsd3 = getGrainSizeDistribution( gss , sizeParam=3 ) # [ Lable, Volume(px3), Size(px), percent passing(%) ]
    gsd4 = getGrainSizeDistribution( gss , sizeParam=4 ) # [ Lable, Volume(px3), Size(px), percent passing(%) ]

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

def relativeBreakage(gsdCurrent,gsdOriginal,maxSize,fracDim):
    gsdOrig, gsdCur, gsdUlt = formatGradationsAndGetUltimate(gsdCurrent,gsdOriginal,maxSize,fracDim)
    Bp = getAreaBetweenGSDs( gsdUlt , gsdOrig )
    B = getAreaBetweenGSDs( gsdCur , gsdOrig )
    Br = B/Bp
    return Br

def getParticleSize(labelledMapForParticleSizeAnalysis):
    numberOfParticles = int(labelledMapForParticleSizeAnalysis.max())

    # Particle size summary columns
    # [0] Index, [1] Volume, [2] Eqsp, [3] Centroidal - max, [4] Centroidal - med, [5] Centroidal - min, [6] and [7] are open
    particleSizeDataSummary = np.zeros( ( numberOfParticles + 1 , 8 ) )
    print( "Starting measurement of particles..." )

    '''
    TODO: This can be parallelized
    '''
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

        caMax = max( feretDims )[ 0 ] # Param 2
        caMin = min( feretDims )[ 0 ] # Param 3
        caMed = statistics.median( feretDims )[ 0 ] # Param 4

        particleSizeDataSummary[particleNum, 3] = caMax
        particleSizeDataSummary[particleNum, 4] = caMed
        particleSizeDataSummary[particleNum, 5] = caMin

    return particleSizeDataSummary  # [ Label, Volume(px3), Size1(px), Size2(px), Size3(px), Size4(px), Size5(ND), Size6(ND)]

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
    print('\tDone')
    return gsdPP # [ Lable, Volume(px), Size(px), percent passing(%) ]

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

def getMorphology(aggregate):
    print("Measuring particle morphology...")
    '''
    Compute some measure of roundness for a particle
    Sphericity - stick to probably ratio of "Feret sizes"
    '''

def measureContactNormalsSpam(aggregate):
    print("\nMeasuring contact normals using SPAM library\n")
    labelledData=aggregate.labelledMap
    binaryData=aggregate.binaryMap
    contactVolume, Z, contactsTable, contactingLabels = slab.labelledContacts(labelledData)

    print("\nMeasuring contact using randomwalker\n")
    ortTabSandRW = slab.contacts.contactOrientationsAssembly(labelledData, binaryData, contactingLabels, watershed="RW")
    tempOrtsRW = np.zeros_like(ortTabSandRW)
    ortOnlySandRW = ortTabSandRW[:,2:5]
    j=0
    for i in range(0,ortTabSandRW.shape[0]):
        if (ortOnlySandRW[i,0]<0):
            ortOnlySandRW[i]=-1*ortOnlySandRW[i]
        if (ortOnlySandRW[i]**2).sum()<=0.999:
            print("Contact deleted - small contact")
        else:
            tempOrtsRW[j] = ortTabSandRW[i]
            j=j+1
    aggregate.contactTableRW = tempOrtsRW[0:j,:]
    np.savetxt("ContactTableRW.csv", aggregate.contactTableRW, delimiter=",")
    print("\nMeasuring contact using watershed\n")
    ortTabSandITK = slab.contacts.contactOrientationsAssembly(labelledData, binaryData, contactingLabels, watershed="ITK")
    tempOrtsITK = np.zeros_like(ortTabSandITK)
    ortOnlySandITK = ortTabSandITK[:,2:5]
    j=0
    for i in range(0,ortTabSandITK.shape[0]):
        if (ortOnlySandITK[i,0]<0):
            ortOnlySandITK[i]=-1*ortOnlySandRW[i]
        elif (ortOnlySandITK[i]**2).sum()<=0.999:
            print("Contact deleted - small contact")
        else:
            tempOrtsITK[j] = ortTabSandITK[i]
            j=j+1
    aggregate.contactTableITK = tempOrtsITK[0:j,:]
    np.savetxt("ContactTableITK.csv", aggregate.contactTableITK, delimiter=",")

def getAreaBetweenGSDs(gsdUp,gsdDown,bins=1000):
    x1 = gsdUp[:,0]
    y1 = gsdUp[:,1]
    x2 = gsdDown[:,0]
    y2 = gsdDown[:,1]

    logX1 = np.log10(x1)
    logX2 = np.log10(x2)

    incrLogX = ( logX1.max() - logX1.min() ) / bins
    logXInterp = np.arange(logX1.min(),logX1.max(),incrLogX)
    y1Interp = np.zeros_like(logXInterp)
    y2Interp = np.zeros_like(logXInterp)
    area = np.zeros_like(logXInterp)

    for i in range(0,bins):
        y1Interp[i] = np.interp(logXInterp[i], logX1, y1 )
        y2Interp[i] = np.interp(logXInterp[i], logX2, y2 )

        if i == 0:
            area[ i ] = ( y1Interp[ i ] - y2Interp[ i ] ) * incrLogX / 2
        elif i == bins - 1:
            area[ i ] = ( y1Interp[ i ] - y2Interp[ i ] ) * incrLogX / 2
        else:
            area[ i ] = ( y1Interp[ i ] - y2Interp[ i ] ) * incrLogX

    totalArea = np.sum(area)

    return totalArea

def formatGradationsAndGetUltimate(gsdCurrent,gsdOriginal,maxSize,fracDim):
    xmax = float(maxSize)
    xmin = float(getStartOfUltimateGradation(xmax,fracDim))
    ymax = float(100)
    ymin = float(0)

    gsdCurrent = gsdCurrent.astype(float)
    gsdOriginal = gsdOriginal.astype(float)

    top = np.array([xmin,ymin]).reshape(1,2)
    bottom = np.array([xmax,ymax]).reshape(1,2)

    top = top.astype( float )
    bottom = bottom.astype( float )

    corGsdCurr = np.append( np.append( top , gsdCurrent , axis = 0 ), bottom, axis=0)
    corGsdOrig = np.append( np.append( top , gsdOriginal, axis = 0 ), bottom, axis=0)
    corGSDUlt = getUltimateGradation(xmax,xmin,fracDim)

    return corGsdOrig, corGsdCurr, corGSDUlt

def getStartOfUltimateGradation(dm,fracDim):
    if VERBOSE: print('\nGetting minimum size of ultimate gradation.')
    d = dm
    small = False

    while small != True:
        pp = (d/dm)**(3-fracDim)
        if VERBOSE: print('\tParticle size is: ' + str(np.round(pp,4)) + 'mm')
        if VERBOSE: print('\tPercentage passing is: ' + str(np.round(pp*100,2)) + '%' )
        if pp > (0.1/100): d*=0.5
        else: small=True

    return xmin

def getUltimateGradation(xmax,xmin,fracDim):
    incr = ( xmax - xmin ) / 100
    x = np.arange(xmin,xmax + incr,incr)
    y = ( ( x / xmax ) ** ( 3 - fracDim ) ) * 100
    x = x.reshape(x.shape[ 0 ] , 1 )
    y = y.rehsape(y.shape[ 0 ] , 1 )
    ultimateGradation = np.append(x,y,axis=1)
    return ultimateGradation.astype(float)

