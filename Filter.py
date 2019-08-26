# -*- coding: utf-8 -*-

"""

"""

# %% Importing libraries

from skimage import filters
from skimage import restoration



#%%

class Filter:

    def __init__(self):
      print("Filter tool activated")
    
    def filterDenoiseNlm(self, aggregate):
      print('Starting filtering')
      inputImage = aggregate.greyLevelMap
      runFilter = True
      patchSize = 3
      patchDistance = 11
      hVal = 0.02
      sigmaEstimate = np.mean(restoration.estimate_sigma(inputImage, multichannel=True))
      while runFilter == True:
        filterImage = restoration.denoise_nl_means(inputImage,
                                                   patch_size=patchSize,
                                                   patch_distance=patchDistance,
                                                   h=hVal,
                                                   multichannel=False,
                                                   fast_mode=True,
                                                   sigma=sigmaEstimate)
        
