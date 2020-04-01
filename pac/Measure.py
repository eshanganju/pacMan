'''
Measure class

'''

import numpy as np
import math
import statistics
import scipy
import spam.label as slab

def gsd(labelledMap):
    gss = getParticleSize( labelledMap )
    gsd1 = getGrainSizeDistribution( gss , sizeParam=1 )
    gsd2 = getGrainSizeDistribution( gss , sizeParam=2 )
    gsd3 = getGrainSizeDistribution( gss , sizeParam=3 )
    return gsd1, gsd2, gsd3

def getParticleSize(labelledMapForParticleSizeAnalysis):
    numberOfParticles = int(labelledMapForParticleSizeAnalysis.max())

    # Particle size summary columns (0) Index, (1) Volume, (2) Eqsp, (3) Centroidal - max, (4) Centroidal - med, (5) Centroidal - min, (6) and (7) are open
    particleSizeDataSummary = np.zeros( ( numberOfParticles + 1 , 8 ) )
    print( "Starting measurement of particles..." )

    '''
    TODO: This can be parallelized
    '''
    for particleNum in range( 1, numberOfParticles + 1 ):
        print( "Computing size of", particleNum, "/", numberOfParticles, "particle" )
        particleSizeDataSummary[particleNum, 0] = particleNum

        # Equivalent sphere diameter
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

        caMax = max( feretDims )[ 0 ]
        caMin = min( feretDims )[ 0 ]
        caMed = statistics.median( feretDims )[ 0 ]

        particleSizeDataSummary[particleNum, 3] = caMax
        particleSizeDataSummary[particleNum, 4] = caMed
        particleSizeDataSummary[particleNum, 5] = caMin

    return particleSizeDataSummary

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
    return gsdPP

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


