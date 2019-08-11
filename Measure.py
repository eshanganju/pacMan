# -*- coding: utf-8 -*-
"""
Outline:
This a part of the PAC code; this modue is for measurement of size, fabric and state parameters of the segmented CT data

---
Features: 
    1. 


---
References: 
    [1]

"""

# %% Importing libraries

# Python libs
import numpy as np
import math
import statistics
import pandas as pd

# Spam libraries
import spam.label as slab

#%% 

class Measure:

    # Initialize
    def __init__(self):
        print("Ruler has been polished: measuring device activated")

    def measureBenchMarkSizeAndNormal(self,aggregate,radii,centres):
        aggregate.benchMarkNumberOfParticles = radii.shape[0]
        aggregate.benchMarkCentres = np.zeros((aggregate.benchMarkNumberOfParticles,3))
        aggregate.benchMarkRadii = np.zeros((aggregate.benchMarkNumberOfParticles,1))
        aggregate.benchMarkCentres = centres
        aggregate.benchMarkRadii = radii

        # Compute grain size distribution
        aggregate.benchMarkGrainSizeDistribution = np.zeros((aggregate.benchMarkNumberOfParticles,2))
        aggregate.benchMarkGrainSizeDistribution[:,0] = aggregate.benchMarkRadii
        aggregate.benchMarkGrainSizeDistribution[:,0] = aggregate.benchMarkGrainSizeDistribution[:,0]*2
        aggregate.benchMarkGrainSizeDistribution = np.sort(aggregate.benchMarkGrainSizeDistribution.view('f8,f8'), order=['f1'], axis=0).view(np.float)
        for i in range(0,aggregate.benchMarkNumberOfParticles):
            aggregate.benchMarkGrainSizeDistribution[i,1]=((sum((aggregate.benchMarkGrainSizeDistribution[0:i+1,0])**3))/sum((aggregate.benchMarkGrainSizeDistribution[:,0])**3))*100

        # Find contact normals of benchMarkData
        aggregate.benchMarkNumberOfContacts = 0
        contactingParticlesOne = []
        contactingParticlesTwo = []
        for i in range(0,aggregate.benchMarkNumberOfParticles-1):
            for j in range(i+1,aggregate.benchMarkNumberOfParticles):
                a = aggregate.benchMarkCentres[i,0]-aggregate.benchMarkCentres[j,0]
                b = aggregate.benchMarkCentres[i,1]-aggregate.benchMarkCentres[j,1]
                c = aggregate.benchMarkCentres[i,2]-aggregate.benchMarkCentres[j,2]
                distanceCentres = (a**2+b**2+c**2)**0.5
                if distanceCentres <= (aggregate.benchMarkRadii[i]+aggregate.benchMarkRadii[j]):
                    aggregate.benchMarkNumberOfContacts = aggregate.benchMarkNumberOfContacts+1
                    contactingParticlesOne.append(i)
                    contactingParticlesTwo.append(j)

        aggregate.benchMarkContactNormal = np.zeros((aggregate.benchMarkNumberOfContacts,5))
        for k in range(0,len(contactingParticlesOne)):
            aggregate.benchMarkContactNormal[k,0] = contactingParticlesOne[k]
            aggregate.benchMarkContactNormal[k,1] = contactingParticlesTwo[k]
            
            a = aggregate.benchMarkCentres[contactingParticlesOne[k],0]-aggregate.benchMarkCentres[contactingParticlesTwo[k],0]
            b = aggregate.benchMarkCentres[contactingParticlesOne[k],1]-aggregate.benchMarkCentres[contactingParticlesTwo[k],1]
            c = aggregate.benchMarkCentres[contactingParticlesOne[k],2]-aggregate.benchMarkCentres[contactingParticlesTwo[k],2]
            d = (a**2+b**2+c**2)**0.5

            # If Z if positive, the sign of all are retained
            if a>=0:
                aggregate.benchMarkContactNormal[k,2] = a/d
                aggregate.benchMarkContactNormal[k,3] = b/d
                aggregate.benchMarkContactNormal[k,4] = c/d
            # If Z is negative, the sign of all are flipped to get positive Z
            elif a<=0:
                aggregate.benchMarkContactNormal[k,2] = -a/d
                aggregate.benchMarkContactNormal[k,3] = -b/d
                aggregate.benchMarkContactNormal[k,4] = -c/d



    def measureParticleSizeDistribution(self,aggregate):

        aggregate.particleSizeDataSummary = np.zeros((aggregate.numberOfParticles+1,6)) # Starting at 1 Index, volume, sphere, Feret min, max and med
    
        print("Starting measurement of particles...")

        for i in range(1,aggregate.numberOfParticles+1):                                # Starting from 1 as 0 is void
            print("Computing size of",i,"/",aggregate.numberOfParticles, "particle")    # Impatience
            # Equivalent sphere diameter
            vol = aggregate.particleList[i].volume
            sphDia = 2*(((3*vol)/(4*math.pi))**(1/3))                                   # Obtain diameter of equivalent sphere
            aggregate.particleList[i].equivalentSphereDiameter = sphDia                 # Assign to particle object

            # Feret diameters
            pointCloud = aggregate.particleList[i].locationData                         # Storing particle point data as a local variable
            aggregate.particleList[i].covarianceMatrix=np.cov(pointCloud.T)             # Covariance matrix (CovMat) of the point cloud
            eigval,eigvec = np.linalg.eig(aggregate.particleList[i].covarianceMatrix)   # computing eigenvalues and igenvecors of the CovMat
            aggregate.particleList[i].eigenValue = eigval                               # Storing eigenvalues as a local variable
            aggregate.particleList[i].eigenVector = eigvec                              # Storing eigenvectors as a local variable
            aggregate.particleList[i].meanZ = np.average(pointCloud[:,0])               # Finding centroid of particle Z
            aggregate.particleList[i].meanY = np.average(pointCloud[:,1])               # Finding centroid of particle Y
            aggregate.particleList[i].meanX = np.average(pointCloud[:,2])               # Finding centroid of particle X
            meanMatrix = np.zeros_like(pointCloud)                                      # Mean matrix of point cloud
            meanMatrix[:,0] = aggregate.particleList[i].meanZ                           # Populating mean matrix of point cloud Z
            meanMatrix[:,1] = aggregate.particleList[i].meanY                           # Populating mean matrix of point cloud Y
            meanMatrix[:,2] = aggregate.particleList[i].meanX                           # Populating mean matrix of point cloud X
            aggregate.particleList[i].centeredLocationData = pointCloud - meanMatrix    # Centering location of centroid over origin
            centeredPointCloud=aggregate.particleList[i].centeredLocationData           # Storing centered point cloud as a local variable
            rotationMatrix=np.zeros((3,3))                                              # Create a rotation matrix
            rotationMatrix[:,0] = eigvec[0]                                             # First column of rotation matrix
            rotationMatrix[:,1] = eigvec[1]                                             # First column of rotation matrix
            rotationMatrix[:,2] = eigvec[2]                                             # First column of rotation matrix
            rotCentPointCloud = (np.matmul(rotationMatrix,centeredPointCloud.T)).T      # Matrix multiplication to get rotated point cloud
            aggregate.particleList[i].rotatedCenteredLocationData = rotCentPointCloud   # Saving rotated, centered particle as an object var
            feretDims = np.zeros((3,1))                                                 # Create local fert dim array
            feretDims[0]=rotCentPointCloud[:,0].max()-rotCentPointCloud[:,0].min()      # Feret dia 1
            feretDims[1]=rotCentPointCloud[:,1].max()-rotCentPointCloud[:,1].min()      # Feret dia 2
            feretDims[2]=rotCentPointCloud[:,0].max()-rotCentPointCloud[:,2].min()      # Feret dia 3
            aggregate.particleList[i].feretMax = max(feretDims)[0]                      # Store largest feret diameter in obect variable
            aggregate.particleList[i].feretMin = min(feretDims)[0]                      # Store smallest feret diameter in obect variable
            aggregate.particleList[i].feretMed = statistics.median(feretDims)[0]        # Store intermediate feret diameter in obect variable

            # Storing data in neat table
            aggregate.particleSizeDataSummary[i][0] = i
            aggregate.particleSizeDataSummary[i][1] = aggregate.particleList[i].volume
            aggregate.particleSizeDataSummary[i][2] = aggregate.particleList[i].equivalentSphereDiameter
            aggregate.particleSizeDataSummary[i][3] = aggregate.particleList[i].feretMax
            aggregate.particleSizeDataSummary[i][4] = aggregate.particleList[i].feretMed 
            aggregate.particleSizeDataSummary[i][5] = aggregate.particleList[i].feretMin 
            # TODO: Add other particle size measure(s)
            # TODO: Use of dataframes for GSD?
            # TODO: Is the eigen vector approach really the feret diameter?

        # Creating grain size distribution list
        tempTable = np.zeros((aggregate.numberOfParticles+1,4))
        tempTable[:,0] = aggregate.particleSizeDataSummary[:,0]
        tempTable[:,1] = aggregate.particleSizeDataSummary[:,1]

        for i in range(2,6):
            tempTable[:,2]= aggregate.particleSizeDataSummary[:,i]
            tempTable = np.sort(tempTable.view('f8,f8,f8,f8'), order=['f2'], axis=0).view(np.float)

            for j in range(1,aggregate.numberOfParticles+1):
                tempTable[j,3]=((sum(tempTable[0:j,1]))/sum(tempTable[:,1]))*100

            if i == 2:
                '''eqSphere'''
                print("eqSphere")
                aggregate.grainSizeDistributionEquivalentSphere = np.zeros((aggregate.numberOfParticles+1,2))
                aggregate.grainSizeDistributionEquivalentSphere[:,0]=tempTable[:,2]
                aggregate.grainSizeDistributionEquivalentSphere[:,1]=tempTable[:,3]

            elif i == 3:
                '''feretMax'''
                print("feretMax")
                aggregate.grainSizeDistributionFeretMax = np.zeros((aggregate.numberOfParticles+1,2))
                aggregate.grainSizeDistributionFeretMax[:,0]=tempTable[:,2]
                aggregate.grainSizeDistributionFeretMax[:,1]=tempTable[:,3]

            elif i== 4:
                '''feretMed'''
                print("feretMed")
                aggregate.grainSizeDistributionFeretMed = np.zeros((aggregate.numberOfParticles+1,2))
                aggregate.grainSizeDistributionFeretMed[:,0]=tempTable[:,2]
                aggregate.grainSizeDistributionFeretMed[:,1]=tempTable[:,3]

            elif i == 5: 
                '''feretMin'''
                print("feretMin")
                aggregate.grainSizeDistributionFeretMin = np.zeros((aggregate.numberOfParticles+1,2))
                aggregate.grainSizeDistributionFeretMin[:,0]=tempTable[:,2]
                aggregate.grainSizeDistributionFeretMin[:,1]=tempTable[:,3]



    def measureMorphology(self,aggregate):
        print("Measuring particle morphology...")
        '''
        Compute some measure of roundness for a particle
        Sphericity - stick to probably ratio of "Feret sizes"
        '''


    def measureContactNormalsSpam(self,aggregate):
        labelledData=aggregate.labelledMap
        binaryData=aggregate.binaryMap
        contactVolume, Z, contactsTable, contactingLabels = slab.labelledContacts(labelledData)             # Uses spam libraries to get contact 
        ortTabSand = slab.contacts.contactOrientationsAssembly(labelledData,binaryData,contactingLabels)    # Uses spam libraries to get contact orientation
        tempOrts = np.zeros_like(ortTabSand)                                                                # Some orientations are bad, cleaning them
        ortOnlySand = ortTabSand[:,2:5]                                                                     # Variable extracts only contacts with enough points
        j=0

        # For some reason the X component of the orientation is kept positive in this code
        # To stick to Z +ve convention of contactOrientation code, we keep Z positive
        for i in range(0,ortTabSand.shape[0]):
            if (ortOnlySand[i,0]<0):
                ortOnlySand[i]=-1*ortOnlySand[i]

            if (ortOnlySand[i]**2).sum()<=0.999:
                print("Contact deleted - small contact")

            else:
                tempOrts[j] = ortTabSand[i]
                j=j+1

        aggregate.contactTable = tempOrts[0:j,:]




