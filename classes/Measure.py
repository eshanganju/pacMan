'''
Measure class

'''

import numpy as np
import math
import statistics
import pandas as pd
import scipy
import spam.label as slab


def measureParticleSize(labelledMapForParticleSizeAnalysis):
    numberOfParticles = int(labelledMapForParticleSizeAnalysis.max())

    # Particle size summary columns (0) Index, (1) Volume, (2) Eqsp, (3) Centroidal - max, (4) Centroidal - med, (5) Centroidal - min, (6) and (7) are open
    particleSizeDataSummary = np.zeros( ( numberOfParticles + 1 , 8 ) )
    print( "Starting measurement of particles..." )

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

        feretDims = np.zeros( ( 3, 1 ) )
        feretDims[ 0 ] = rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 0 ].min()
        feretDims[ 1 ] = rotCentPointCloud[ :, 1 ].max() - rotCentPointCloud[ :, 1 ].min()
        feretDims[ 2 ] = rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 2 ].min()

        feretMax = max( feretDims )[ 0 ]
        feretMin = min( feretDims )[ 0 ]
        feretMed = statistics.median( feretDims )[ 0 ]

        particleSizeDataSummary[particleNum, 3] = feretMax
        particleSizeDataSummary[particleNum, 4] = feretMed
        particleSizeDataSummary[particleNum, 5] = feretMin

    return particleSizeDataSummary

def computeGrainSizeDistribution(particleSizeSummary):
    outputDF = particleSizeSummary.copy()
    totalVol = outputDF['vol'].sum()
    outputDF['mass%'] = outputDF['vol']/totalVol*100

    # EQSP
    outputDF.sort_values('eqsp',inplace=True)
    outputDF.reset_index(drop=True)
    outputDF['eqspPP%'] = outputDF['mass%'].cumsum()

    # caMax
    outputDF.sort_values('caMax',inplace=True)
    outputDF.reset_index(drop=True)
    outputDF['caMaxPP%'] = outputDF['mass%'].cumsum()

    # caMed
    outputDF.sort_values('caMed',inplace=True)
    outputDF.reset_index(drop=True)
    outputDF['caMedPP%'] = outputDF['mass%'].cumsum()

    # caMin
    outputDF.sort_values('caMin',inplace=True)
    outputDF.reset_index(drop=True)
    outputDF['caMinPP%'] = outputDF['mass%'].cumsum()

    return outputDF

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

# Morphology
def measureMorphology(aggregate):
    print("Measuring particle morphology...")
    '''
    Compute some measure of roundness for a particle
    Sphericity - stick to probably ratio of "Feret sizes"
    '''

# REV size analysis
def revSizeAnalysis( aggregate ):
    print('\nStarting REV size analysis')
    '''
    Check different rev sizes and get GSD and other stuff
    '''

# Contact Normals
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


