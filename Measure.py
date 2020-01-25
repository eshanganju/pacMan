# -*- coding: utf-8 -*-
"""
Outline:
This a part of the PAC code; this module is for measurement of size, fabric and state parameters of the segmented CT data

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
#
#import spam.label as slab

#%%  Measure class and methods CalmPatientFocused

class Measure:

    # Initialize
    def __init__(self):
        print("Measurer activated")

    def measureBenchMarkSizeAndNormal(self,aggregate,radii,centres):      
        aggregate.benchMarkNumberOfParticles = radii.shape[0]        
        aggregate.benchMarkCentres = np.zeros((aggregate.benchMarkNumberOfParticles,3))       
        aggregate.benchMarkRadii = np.zeros((aggregate.benchMarkNumberOfParticles,1))        
        aggregate.benchMarkCentres = centres       
        aggregate.benchMarkRadii = radii

        # Compute grain size distribution
        aggregate.benchMarkGrainSizeDistribution = np.zeros((aggregate.benchMarkNumberOfParticles,2))        
        aggregate.benchMarkGrainSizeDistribution[:,0] = (aggregate.benchMarkRadii)*2        
        aggregate.benchMarkGrainSizeDistribution = np.sort(aggregate.benchMarkGrainSizeDistribution.view('f8,f8'), order=['f1'],axis=0).view(np.float)
        
        for i in range(0,aggregate.benchMarkNumberOfParticles):        
            a = sum((aggregate.benchMarkGrainSizeDistribution[0:i+1,0])**3)           
            b = sum((aggregate.benchMarkGrainSizeDistribution[:,0])**3)            
            aggregate.benchMarkGrainSizeDistribution[i,1]=(a/b)*100

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
              
        np.savetxt("benchMarkContactNormals.csv",aggregate.benchMarkContactNormal,delimiter=",")      
        np.savetxt("benchMarkGSD.csv",aggregate.benchMarkGrainSizeDistribution,delimiter=",")
                

    # Particle Size distribution
    def measureParticleSizeDistribution( self, aggregate ):
        aggregate.particleSizeDataSummary = np.zeros( ( aggregate.numberOfParticles + 1 , 6 ) )
        print( "Starting measurement of particles..." )

        for i in range( 1, aggregate.numberOfParticles + 1 ):                                          
            print( "Computing size of", i, "/", aggregate.numberOfParticles, "particle" )         
            # Equivalent sphere diameter
            vol = aggregate.particleList[ i ].volume        
            sphDia = 2 * ( ( ( 3 * vol ) / ( 4 * math.pi ) ) ** ( 1 / 3 ) )                                       
            aggregate.particleList[ i ].equivalentSphereDiameter = sphDia            


            # Feret diameters
            pointCloud = aggregate.particleList[ i ].locationData                   
            aggregate.particleList[ i ].covarianceMatrix = np.cov( pointCloud.T )      
            eigval, eigvec = np.linalg.eig( aggregate.particleList[ i ].covarianceMatrix )

            aggregate.particleList[ i ].eigenValue = eigval                            
            aggregate.particleList[ i ].eigenVector = eigvec                          
            
            aggregate.particleList[ i ].meanZ = np.average( pointCloud[ :, 0 ] )           
            aggregate.particleList[ i ].meanY = np.average( pointCloud[ :, 1 ] )              
            aggregate.particleList[ i ].meanX = np.average( pointCloud[ :, 2 ] )             

            meanMatrix = np.zeros_like( pointCloud )                                   

            meanMatrix[ :, 0 ] = aggregate.particleList[ i ].meanZ                       
            meanMatrix[ :, 1 ] = aggregate.particleList[ i ].meanY                      
            meanMatrix[ :, 2 ] = aggregate.particleList[ i ].meanX                     

            aggregate.particleList[ i ].centeredLocationData = pointCloud - meanMatrix 

            centeredPointCloud = aggregate.particleList[ i ].centeredLocationData        
            
            rotationMatrix = np.zeros( ( 3, 3 ) )                                          

            rotationMatrix[ :, 0 ] = eigvec[ 0 ]                                        
            rotationMatrix[ :, 1 ] = eigvec[ 1 ]                                       
            rotationMatrix[ :, 2 ] = eigvec[ 2 ]                                            

            rotCentPointCloud = ( np.matmul( rotationMatrix, centeredPointCloud.T ) ).T     
            aggregate.particleList[ i ].rotatedCenteredLocationData = rotCentPointCloud  

            feretDims = np.zeros( ( 3, 1 ) )                                                
            feretDims[ 0 ] = rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 0 ].min()     
            feretDims[ 1 ] = rotCentPointCloud[ :, 1 ].max() - rotCentPointCloud[ :, 1 ].min()     
            feretDims[ 2 ] = rotCentPointCloud[ :, 0 ].max() - rotCentPointCloud[ :, 2 ].min()     

            aggregate.particleList[ i ].feretMax = max( feretDims )[ 0 ]                     
            aggregate.particleList[ i ].feretMin = min( feretDims )[ 0 ]                     
            aggregate.particleList[ i ].feretMed = statistics.median( feretDims )[ 0 ]       


            # Storing data in neat table
            aggregate.particleSizeDataSummary[ i ][ 0 ] = i            
            aggregate.particleSizeDataSummary[ i ][ 1 ] = aggregate.particleList[ i ].volume            
            aggregate.particleSizeDataSummary[ i ][ 2 ] = aggregate.particleList[ i ].equivalentSphereDiameter            
            aggregate.particleSizeDataSummary[ i ][ 3 ] = aggregate.particleList[ i ].feretMax            
            aggregate.particleSizeDataSummary[ i ][ 4 ] = aggregate.particleList[ i ].feretMed             
            aggregate.particleSizeDataSummary[ i ][ 5 ] = aggregate.particleList[ i ].feretMin 
            

        np.savetxt( "DataGSD.csv", aggregate.particleSizeDataSummary, delimiter = "," )        


    # Morphology
    def measureMorphology(self,aggregate):
        print("Measuring particle morphology...")
        '''
        Compute some measure of roundness for a particle
        Sphericity - stick to probably ratio of "Feret sizes"
        '''

    # Contact Normals
#    def measureContactNormalsSpam(self,aggregate):
#        print("\nMeasuring contact normals using SPAM library\n")
#        labelledData=aggregate.labelledMap
#        binaryData=aggregate.binaryMap
#        
#        contactVolume, Z, contactsTable, contactingLabels = slab.labelledContacts(labelledData)
#
#        # Random walker - contact
#        print("\nMeasuring contact using randomwalker\n")
#        ortTabSandRW = slab.contacts.contactOrientationsAssembly(labelledData,
#                                                                 binaryData,
#                                                                 contactingLabels,
#                                                                 watershed="RW")
#        tempOrtsRW = np.zeros_like(ortTabSandRW)
#        ortOnlySandRW = ortTabSandRW[:,2:5]
#        j=0
#        # For some reason the X component of the orientation is kept positive in this code
#        # To stick to Z +ve convention of contactOrientation code, we keep Z positive
#        for i in range(0,ortTabSandRW.shape[0]):
#            if (ortOnlySandRW[i,0]<0):
#                ortOnlySandRW[i]=-1*ortOnlySandRW[i]
#
#            if (ortOnlySandRW[i]**2).sum()<=0.999:
#                print("Contact deleted - small contact")
#
#            else:
#                tempOrtsRW[j] = ortTabSandRW[i]
#                j=j+1
#        aggregate.contactTableRW = tempOrtsRW[0:j,:]
#        np.savetxt("ContactTableRW.csv", aggregate.contactTableRW, delimiter=",")
#
#        # Watershed - contact
#        print("\nMeasuring contact using watershed\n")
#        ortTabSandITK = slab.contacts.contactOrientationsAssembly(labelledData,
#                                                                  binaryData,
#                                                                  contactingLabels,
#                                                                  watershed="ITK")
#        tempOrtsITK = np.zeros_like(ortTabSandITK)
#        ortOnlySandITK = ortTabSandITK[:,2:5]
#        j=0
#        # For some reason the X component of the orientation is kept positive in this code
#        # To stick to Z +ve convention of contactOrientation code, we keep Z positive
#        for i in range(0,ortTabSandITK.shape[0]):
#            if (ortOnlySandITK[i,0]<0):
#                ortOnlySandITK[i]=-1*ortOnlySandRW[i]
#
#            if (ortOnlySandITK[i]**2).sum()<=0.999:
#                print("Contact deleted - small contact")
#
#            else:
#                tempOrtsITK[j] = ortTabSandITK[i]
#                j=j+1
#        aggregate.contactTableITK = tempOrtsITK[0:j,:]
#        np.savetxt("ContactTableITK.csv", aggregate.contactTableITK, delimiter=",")
#
#    def measureContactNormalGeneral(self,aggregate):
#      print("Measuring contact normals - non SPAM - using pixel data")
#      """
#      1. Contact determination for pixel data
#      2. Contact normal from PCA - probabilistic distrbution
#      """
#
#    def measureContactNormalLevelSet(self,aggregate):
#      print("Measuring contact normals - using level set data")

