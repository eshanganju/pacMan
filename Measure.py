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

# Spam libraries
import spam.label as slab

#%% 

class Measure:

    # Initialize
    def __init__(self):
        print("Ruler has been polished, measuring device activated")


    def measureParticleSizeDistribution(self,aggregate):
        '''
        TODO:
        Create a pandas dataframe with diffrent sizes as column headings
        populate that as you compute the particle size parameter
        Advantage would be easier to see and maybe better functionality?
        Plus you get to learn the damn thing
        '''
        print("Starting")
        for i in range(1,aggregate.numberOfParticles+1):                                # Starting from 1 as 0 is void

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
            print("Computing size of",i,"/",aggregate.numberOfParticles, "particle")    # Impatience

            # TODO: Add other particle size measure(s)
            # TODO: Is the eigen vector approach really the feret diameter?
        '''
        Work with dataFrame to get GSD
        '''



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
        for i in range(0,ortTabSand.shape[0]):
            if (ortOnlySand[i]**2).sum()<=0.999:
                print("Contact deleted - small contact")
            else:
                tempOrts[j] = ortTabSand[i]
                j=j+1
        aggregate.contactTable = tempOrts[0:j,:]




