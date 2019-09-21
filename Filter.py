# -*- coding: utf-8 -*-

"""
This is the filtering class - contains filter tools. 
2019-08-27: Coded Non-local Means filter loop - single loop, no input from users

"""

# Importing libraries

from skimage import restoration
import numpy as np
import matplotlib.pyplot as plt
import time

class Filter:

    def __init__(self):
      print("Filter activated")
    
    def filterDenoiseNlm(self, aggregate):
     
      # Adding local variables
      inputImage = aggregate.greyLevelMap
      filteredImage = np.zeros_like(inputImage)
      noiseRemoved = np.zeros_like(filteredImage)
      numHistPts = inputImage.shape[0]*inputImage.shape[1]*inputImage.shape[2]
      his_x = np.arange(0.0005,1,0.001)
      runFilter = True
      
      # Inital settings chosen based on parameteric study 
      patchSize = 5
      patchDistance = 11
      hVal = 0.05
      
      # The SD of the input noisy image
      sigmaEstimate = np.mean(restoration.estimate_sigma(inputImage, multichannel=True))
      
      print('\n\nIntial parameters for filter----------------------*\n')
      print('Patch size: ', patchSize)
      print('Patch distance: ', patchDistance)
      print('Cut-off pixel intensity: ', hVal)

      while runFilter == True:
        
        print('\n\nStarting filtering loop ------------------------*')
        start_time = time.time()

        filteredImage = restoration.denoise_nl_means(inputImage, 
                                                     patch_size=patchSize, 
                                                     patch_distance=patchDistance,
                                                     h=hVal, 
                                                     multichannel=False, 
                                                     fast_mode=True, 
                                                     sigma=sigmaEstimate)
                                                     
        noiseRemoved = abs(inputImage-filteredImage) 
        
        listNoisy = inputImage.reshape((numHistPts,1))
        listClean = filteredImage.reshape((numHistPts,1))   
        
        histNoisy = np.histogram(listNoisy, bins=1000, range=(0,1)) 
        histClean = np.histogram(listClean, bins=1000, range=(0,1))
        
        ax1 = plt.subplot(221)
        ax1.imshow(inputImage[inputImage.shape[0]//2],cmap='Greys_r')
        
        ax2 = plt.subplot(222)
        ax2.imshow(filteredImage[filteredImage.shape[0]//2],cmap='Greys_r')
        
        ax3 = plt.subplot(223)
        ax3.imshow(noiseRemoved[noiseRemoved.shape[0]//2],cmap='Greys_r')
        
        ax4 = plt.subplot(224)
        ax4.plot(his_x,histNoisy[0])
        ax4.plot(his_x,histClean[0])
        
        timeTakenThisLoop = (time.time() - start_time)
        print("\n--- %s seconds ---" %timeTakenThisLoop)
        plt.show()
       
        # Check user: 
        answer = input("\n\nIs the quality of filtered image good [y/n]: ")
        
        if answer == 'n':
          print("Enter new parameters: ")
          patchSize = int(input("Patch size (integer 3-7): "))
          patchDistance = int(input("Patch distance (integer 10-15): "))
          hVal = float(input("Cut-off pixel value (float 0-1): "))
          runFilter = True

        else:
          print("Done! Filtered image quality acceptable----------*")
          runFilter = False
      
      aggregate.filteredGreyLevleMap = filteredImage

      return patchSize, patchDistance, hVal

