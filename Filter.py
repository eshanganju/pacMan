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
    
    def filterDenoiseNlm(self, aggregate):

        # Local variables for central CS of maps
        
        inputImage = aggregate.greyLevelMap[ aggregate.greyLevelMap.shape[ 0 ] // 2 ]
        
        filteredImage = np.zeros_like( inputImage )
        
        noiseRemoved = np.zeros_like( filteredImage )

        # Local variables for the entire map 
        
        inputMap = aggregate.greyLevelMap
        
        filteredMap = np.zeros_like( inputMap )

        # Histograms
        
        numHistPts = inputImage.shape[0]*inputImage.shape[1]
        
        numHistPtsMap = inputMap.shape[0]*inputMap.shape[1]*inputMap.shape[2]
        
        his_x = np.arange( 0, 65536, 1 )

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
        
        runFilter = True
        
        while runFilter == True:
            
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
            
            histNoisy = np.histogram(listNoisy, bins = 65536, range = ( 0, 65535 ) ) 
            
            histClean = np.histogram(listClean, bins = 65536, range = ( 0, 65535 ) )
            
            histNoisy = histNoisy / ( histNoisy[0].sum() ) * 100
            
            histClean = histClean / ( histClean[0].sum() ) * 100

            # Figures
            plt.figure()
            plt.imshow(noiseRemoved, cmap= 'Greys_r' )
            plt.draw() # draw the plot
            plt.pause( 10 ) # show it for 30 seconds
            plt.close()
            
            # Figures
            plt.figure()
            plt.imshow(filteredImage, cmap= 'Greys_r' )
            plt.draw() # draw the plot
            plt.pause( 10 ) # show it for 30 seconds
            plt.close()

            # Figures - Histogram
            plt.figure()
            plt.plot(his_x, histNoisy[ 0 ])
            plt.plot(his_x, histClean[ 0 ])
            plt.grid()
            plt.draw() # draw the plot
            plt.pause( 10 ) # show it for 30 seconds
            plt.close()
        
            timeTakenThisLoop = ( time.time() - start_time )
            
            print( "\n--- %s seconds ---" %round( timeTakenThisLoop ) )

            answer = input( "\n\nAre filter parameters suitable (y/n)?:" )
            
            if answer == 'y':
                runFilter = False
            
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
        histNoisy = np.histogram( listNoisy, bins = 65536, range = ( 0, 65535 ) ) 
        
        histClean = np.histogram( listClean, bins = 65536, range = ( 0, 65535 ) )
       
        histNoisy = histNoisy / ( histNoisy[ 0 ].sum() ) * 100
       
        histClean = histClean / ( histClean[ 0 ].sum() ) * 100
    
    
        # Figures - Histogram
        plt.figure()
        
        plt.plot(his_x, histNoisy[ 0 ])
        
        plt.plot(his_x, histClean[ 0 ])
        
        plt.grid()
        
        plt.draw() # draw the plot
        
        plt.pause( 30 ) # show it for 30 seconds
        
        plt.close()
        
        print( '.\n.\n.\nFilter for entire grey level map completed--------*' )
        
        timeTakenThisLoop = ( time.time() - start_time )
        
        print( "\n--- Time taken: %s seconds ---\n" % round( timeTakenThisLoop ) )

        # Update aggregate object
        aggregate.filteredGreyLevelMap = filteredMap
        
        aggregate.imageNoise = aggregate.greyLevelMap - aggregate.filteredGreyLevelMap

        # Keeping records
        tiffy.imsave('Filtered.tiff', aggregate.filteredGreyLevelMap)
        f = open("NonLocalMeansFilterParameters.txt","w+")
        f.write("NLM filter parameters---------------------*\n") 
        f.write("Patch size = %f\n" % patchSize) 
        f.write("Patch distance = %f\n" % patchDistance)
        f.write("Cut-off intensity = %f\n" % hVal)
        f.write("Estimated sigma - CS = %f\n" % sigmaEstimateImage)
        f.write("Estimated sigma - Map = %f\n" % sigmaEstimateMap)
        f.close()

