Approach
  Make particle data using kalipshera
  Add partial volume effect
  Add gaussian blurring
  Add random noise
  ---
  Find effective ways to clean noise
  Check for best performing settings
  Code Filter.py class with the best methods

Knead to know
  The median filter takes a 2D/3D patch over the image and finds median gl of the patch
    The gl at center of the patch is replaced by the median value
  The non-local-means filter exploits self similatity of images [2]
    The non-local means algorithm replaces the value of a pixel by an average of a selection of 
    other pixels values: small patches centered on the other pixels are compared to the patch 
    centered on the pixel of interest, and the average is performed only for pixels that have 
    patches close to the current patch. As a result, this algorithm can restore well textures, 
    that would be blurred by other denoising algorithm.
        
Issues: 
  The nlm filters make the gl peaks narrower but shift the peaks for the particles to smaller gl
  Need to find a way of checking accuracy
    Different in baseline and non-filtered image 
      Mean of the difference of image with PVE and after filters is 0.001395 
    
Conclusions: 
  The non local mean filter does the best job for cleaning the data - the effect of patch
  The  median filter is also good but more noise on the void part
  Median filter does a goood job of identifying the peak locaiton for thye particles
  The difference between the images - original and after sharpening is shown below and saved.

References: 
[1] Vlahinić, Ivan, Edward Andò, Gioacchino Viggiani, and José E. Andrade. 2014. “Towards a More Accurate Characterization of Granular Media: Extracting Quantitative Descriptors from Tomographic Images.” Granular Matter 16 (1): 9–21. https://doi.org/10.1007/s10035-013-0460-6.
[2] Buades, Antoni, Bartomeu Coll, and Jean Michel Morel. 2008. “Nonlocal Image and Movie Denoising.” International Journal of Computer Vision 76 (2): 123–39. https://doi.org/10.1007/s11263-007-0052-1.

Links: 
[1] NLM: https://www.youtube.com/watch?v=9tUns4HYtcw
[2] NLM math: https://www.youtube.com/watch?v=k8hS6uTz-Uc
