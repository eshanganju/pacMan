# -*- coding: utf-8 -*-

"""
This is the filtering class - contains filter tools. 
2019-08-27: Coded Non-local Means filter loop - single loop, no input from users

"""

# Importing libraries

from skimage import restoration
import skimage.external.tifffile as tiffy
import numpy as np
import matplotlib.pyplot as plt
import time

class Filter:

    def __init__(self):
        
        print("Filter activated")
        
    def checkForFiltration(self, aggregate):
        
        checkForFiltrationData = input( 'Is filtered data available? [y/n]: ' )
        
        if checkForFiltrationData == 'y':
            
            print('Default data locaiton is: C:/Users/eganj/Google Drive/(02) Research/(07) EIDPS/(02) Data/(05) Tomo/(02) Physics - nCT/(05) 1D compression study/(02) REV analysis/density based binarization/OTC/superCube/Filtered.tiff')
            defaultDataLocationCorrect = input( 'is this correct? [y/n]: ' )
            
            if defaultDataLocationCorrect == 'y':
                print('Reading GLI data from default location')
                aggregate.filteredGreyLevelMap = tiffy.imread( 'C:/Users/eganj/Google Drive/(02) Research/(07) EIDPS/(02) Data/(05) Tomo/(02) Physics - nCT/(05) 1D compression study/(02) REV analysis/density based binarization/OTC/superCube/Filtered.tiff' )
                # Add a plotting stage
                
            else:
                newDataLocation = input( 'Enter filtered tiff stack location: ' )
                print('Reading GLI data from new location')
                aggregate.filteredGreyLevelMap = tiffy.imread(newDataLocation)
                
            aggregate.imageNoise = aggregate.greyLevelMap - aggregate.filteredGreyLevelMap
            
            # Keeping records
            filteredFileName = aggregate.fileName + '-filtered.tiff'
            noiseFileName =  aggregate.fileName + '-noise.tiff'
            tiffy.imsave(filteredFileName, aggregate.filteredGreyLevelMap)
            tiffy.imsave(noiseFileName, aggregate.imageNoise)
            
        else:
            aggregate.filteredGreyLevelMap = self.filterDenoiseNlm(aggregate.fileName, aggregate.greyLevelMap, aggregate.GLIMax)
            aggregate.imageNoise = aggregate.greyLevelMap - aggregate.filteredGreyLevelMap
            
            # Keeping records
            filteredFileName = aggregate.fileName + '-filtered.tiff'
            noiseFileName =  aggregate.fileName + '-noise.tiff'
            tiffy.imsave(filteredFileName, aggregate.filteredGreyLevelMap)
            tiffy.imsave(noiseFileName, aggregate.imageNoise)

            
    def filterDenoiseNlm(self, sandName, gli, gliMax):

        # Local variables for central CS of maps
        inputImage = gli[ gli.shape[ 0 ] // 2 ] 
        filteredImage = np.zeros_like( inputImage )        
        noiseRemoved = np.zeros_like( filteredImage )

        
        # Local variables for the entire map        
        inputMap = gli       
        filteredMap = np.zeros_like( inputMap )

       
        # Histograms      
        numHistPts = inputImage.shape[0]*inputImage.shape[1]       
        numHistPtsMap = inputMap.shape[0]*inputMap.shape[1]*inputMap.shape[2]      
        his_x = np.arange( 0, gliMax + 1, 1 )

      
        # The SD of the input noisy image      
        sigmaEstimateImage = np.mean(restoration.estimate_sigma(inputImage, multichannel=True))      
        sigmaEstimateMap = np.mean(restoration.estimate_sigma(inputMap, multichannel=True))
             
        
        # Inital settings chosen based on parameteric study         
        patchSize = 3
        patchDistance = 5        
        hVal = 450

        print('\n\nIntial parameters for filter----------------------*\n')
        print('Patch size: ', patchSize)
        print('Patch distance: ', patchDistance)
        print('Cut-off pixel intensity: ', round(hVal))

        
        # Loop for checking filtering parameters; run on CS of data        
        runSliceFilter = True        
        while runSliceFilter == True:           
            print('\n\nNew filtering loop started on cross-section--------*\n\n\n')            
            start_time = time.time()            
            filteredImage = restoration.denoise_nl_means(inputImage,
                                                         patch_size = patchSize, 
                                                         patch_distance = patchDistance, 
                                                         h = hVal, 
                                                         multichannel = False, 
                                                         fast_mode = True, 
                                                         sigma = sigmaEstimateImage)
            
            noiseRemoved = abs( inputImage - filteredImage )             
            listNoisy = inputImage.reshape( ( numHistPts, 1 ) )           
            listClean = filteredImage.reshape( ( numHistPts, 1 ) )

            # Histogram            
            histNoisy = np.histogram(listNoisy, bins = GLIMax + 1, range = ( 0, GLIMax ) ) 
            histClean = np.histogram(listClean, bins = GLIMax + 1, range = ( 0, GLIMax ) )            
            histNoisy = histNoisy / ( histNoisy[0].sum() ) * 100            
            histClean = histClean / ( histClean[0].sum() ) * 100

            # Figures
            plt.figure()
            plt.imshow(noiseRemoved, cmap= 'Greys_r' )
            plt.draw() # draw the plot
            plt.savefig('NoiseRemoved-CenterSlice.png')
            plt.pause( 5 ) # show it for 5 seconds
            plt.close()
            
            # Figures
            plt.figure()
            plt.imshow(filteredImage, cmap= 'Greys_r' )
            plt.draw() # draw the plot
            plt.savefig('filteredImage-CenterSlice.png')
            plt.pause( 5 ) # show it for 5 seconds
            plt.close()

            # Figures - Histogram
            plt.figure()
            plt.plot(his_x, histNoisy[ 0 ])
            plt.plot(his_x, histClean[ 0 ])
            plt.grid()
            plt.draw() # draw the plot
            plt.savefig('histogram-CenterSlice.png')
            plt.pause( 5 ) # show it for 5 seconds
            plt.close()
        
            timeTakenThisLoop = ( time.time() - start_time )            
            print( "\n--- %s seconds ---" %round( timeTakenThisLoop ) )
            answer = input( "\n\nCheck files - are filter parameters suitable (y/n)?:" )
            
            if answer == 'y':
                runSliceFilter = False
            
            else:
                print( "Enter new parameters: \n\n" )                
                patchSize = int( input( "Patch size (pixel): " ) )               
                patchDistance = int( input( "Patch distance (pixel):" ) )               
                hVal = float( input( "Cut-off pixel intensity: " ) )

        
        # Filtering entire grey level map with chosen parameters
        print( '\n\nFilter for entire grey level map started--------*' )
        start_time = time.time()        
        filteredMap = restoration.denoise_nl_means(inputMap, 
                                                   patch_size = patchSize, 
                                                   patch_distance = patchDistance,
                                                   h = hVal, 
                                                   multichannel = False, 
                                                   fast_mode = True, 
                                                   sigma = sigmaEstimateMap)
        
        listNoisy = inputMap.reshape( ( numHistPtsMap, 1 ) )        
        listClean = filteredMap.reshape( ( numHistPtsMap, 1 ) )

        
        # Histogram
        histNoisy = np.histogram( listNoisy, bins = gliMax + 1, range = ( 0, gliMax ) )         
        histClean = np.histogram( listClean, bins = gliMax + 1, range = ( 0, gliMax ) )       
        histNoisy = histNoisy / ( histNoisy[ 0 ].sum() ) * 100      
        histClean = histClean / ( histClean[ 0 ].sum() ) * 100
    
    
        # Figures - Histogram
        plt.figure()        
        plt.plot(his_x, histNoisy[ 0 ])        
        plt.plot(his_x, histClean[ 0 ])        
        plt.grid()
        plt.draw() # draw the plot 
        volumeHistogramName = sandName+'GLIhistogram-volume.png'
        plt.savefig(volumeHistogramName)
        plt.pause( 5 ) # show it for 5 seconds       
        plt.close()
        
        print( '.\n.\n.\nFilter for entire grey level map completed--------*' )       
        timeTakenThisLoop = ( time.time() - start_time )        
        print( "\n--- Time taken: %s seconds ---\n" % round( timeTakenThisLoop ) )
        
        # Updating filtered parameters
        f = open("NonLocalMeansFilterParameters.txt","w+")
        f.write("NLM filter parameters---------------------*\n") 
        f.write("Patch size = %f\n" % patchSize) 
        f.write("Patch distance = %f\n" % patchDistance)
        f.write("Cut-off intensity = %f\n" % hVal)
        f.write("Estimated sigma - CS = %f\n" % sigmaEstimateImage)
        f.write("Estimated sigma - Map = %f\n" % sigmaEstimateMap)
        f.write("Time taken for filter (s) = %f\n" % round(timeTakenThisLoop))
        f.close()

        return filteredMap


