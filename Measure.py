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
from sklearn.decomposition import PCA

# Spam libraries
import spam.label as slab

#%% 

class Measure:

    # Initialize
    def __init__(self):
        print("Ruler has been polished, measuring device activated")


    def measureParticleSizeDistribution(self,aggregate):

        for i in range(1,aggregate.numberOfParticles+1):
            vol = aggregate.particleList[i].volume
            sphDia = 2*(((3*vol)/(4*math.pi))**(1/3))                                   # Obtain diameter of equivalent sphere
            aggregate.particleList[i].equivalentSphereDiameter = sphDia                 # Assign to particle object

            # Finding Feret diameters
            pointCloud = aggregate.particleList[i].locationData                         # Storing particle point data as a local variable
            aggregate.particleList[i].covarianceMatrix=np.cov(pointCloud.T)             # 
            eigval,eigvec = np.linalg.eig(aggregate.particleList[i].covarianceMatrix)
            aggregate.particleList[i].eigenValue = eigval
            aggregate.particleList[i].eigenVector = eigvec
            aggregate.particleList[i].meanZ = np.average(pointCloud[:,0])
            aggregate.particleList[i].meanY = np.average(pointCloud[:,1])
            aggregate.particleList[i].meanX = np.average(pointCloud[:,2])
            meanMatrix = np.zeros_like(pointCloud)
            meanMatrix[:,0] = aggregate.particleList[i].meanZ
            meanMatrix[:,1] = aggregate.particleList[i].meanY
            meanMatrix[:,2] = aggregate.particleList[i].meanX
            aggregate.particleList[i].centeredLocationData = pointCloud - meanMatrix
            centPointCloud=aggregate.particleList[i].centeredLocationData
            '''
            Rotate centPointCloud from currrent basis to the eigen vector basis
            Get max and min values of points on each axes
            Differenc of max and min values along each axes are the three feret diameters
            '''

    def measureMorphology(self,aggregate):
        '''
        Compute some measure of roundness for a particle
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




