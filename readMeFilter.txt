Objective:
  Flter noise from image to have more well defined gl intensity distributions. 
  Check for fidelity with the original image with PVE only 

Approach:
  Using Kalisphera: 
    Get particle data 
    Add partial volume effect
    Add gaussian blurring
    Add random noise
  ---
  Find effective ways to clean noise using different filters
  Check for best performing settings
  Code Filter.py class with the best methods
  The trial filter files are named as alphaFilter.py and betaFilter.py, utill the final class
    is named as Filter.py 
  
Knead to know
  In literature [1], the main methods listed are [gaussian], [median] and [non-local means]
    filter to clean the noise in the image. The Gaussian filter does a local weighted
    averaging of the image over a patch, using gaussian distribution as the weights.
    This burs out the image - no good.
  The median filter takes a 2D/3D patch over the image and finds median gl of the patch
    The gl at center of the patch is replaced by the median value. This tends to preserve
    the edge of the particles better (depending on the size of patch chosen). 
  The non-local-means filter exploits self similatity of images [2]
    The non-local means algorithm replaces the value of a pixel by a weighted average of 
    a selection of other pixels values: small patches centered on the other 
    pixels are compared to the patch centered on the pixel of interest, and the 
    average is performed only for pixels that have patches close to the current 
    patch. As a result, this algorithm can restore well textures, that would be 
    blurred by other denoising algorithm.
        
Findings: 
  1. nlm filters make the gl peaks narrower but shift left the peaks for the 
  2. Need to find a way of checking accuracy - images look better but needs to
    Different in baseline and non-filtered image 
    Mean of the difference of image with PVE and after filters is 0.001395 
    
Conclusions: 
  The non local mean filter does the best job for cleaning the data - the effect of patch
  The  median filter is also good but more noise on the void part
  Median filter does a goood job of identifying the peak locaiton for thye particles
  The difference between the images is saved 

References: 
[1] Vlahinić, Ivan, Edward Andò, Gioacchino Viggiani, and José E. Andrade. 2014. “Towards a More Accurate Characterization of Granular Media: Extracting Quantitative Descriptors from Tomographic Images.” Granular Matter 16 (1): 9–21. https://doi.org/10.1007/s10035-013-0460-6.
[2] Buades, Antoni, Bartomeu Coll, and Jean Michel Morel. 2008. “Nonlocal Image and Movie Denoising.” International Journal of Computer Vision 76 (2): 123–39. https://doi.org/10.1007/s11263-007-0052-1.

Links: 
[1] NLM: https://www.youtube.com/watch?v=9tUns4HYtcw
[2] NLM math: https://www.youtube.com/watch?v=k8hS6uTz-Uc

