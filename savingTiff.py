"""
Making my stupid ass figure out how to store numpy stack as atiff image
"""
# General [Tso]
import numpy as np
import matplotlib.pyplot as plt
import skimage.external.tifffile as tiffy
import math
import scipy

# Spam
import spam.datasets as sdata
import spam.kalisphera as skali

pixelSize = 30.e-6                                                      
blurSTD = 0.8                                                          
noiseSTD = 0.03                                                       
boxSizeDEM, centres, radii = sdata.loadDEMboxsizeCentreRadius()      
rMax = np.amax(radii)                                               
boxSize = boxSizeDEM + 3 * rMax                                         
centres[:, :] = centres[:, :] + 1.5 * rMax                             
boxSize = int(math.ceil(np.max(boxSize[:]) / pixelSize))
centres = centres / pixelSize
radii = radii / pixelSize
Box = np.zeros((boxSize, boxSize, boxSize), dtype="<f8")                
skali.makeSphere(Box, centres, radii)                                   
Box1 = Box
Box1[np.where(Box1 > 1.0)] = 1.0                                        
Box1[np.where(Box1 < 0.0)] = 0.0                                        
Box2 = Box1 * 0.5                                                       
Box2 = Box2 + 0.25                                                      
Box3 = scipy.ndimage.filters.gaussian_filter(Box2, sigma=blurSTD)       
Box4 = np.random.normal(Box3, scale=noiseSTD)                           

Box4A = Box4

# Save file as imageJ stack
tiffy.imsave('Box4A.tiff',Box4A)

# Load the file into var and show
Box4B = tiffy.imread('Box4A.tiff')

# Plot the written and read data to compare
plt.subplot(131)
plt.imshow(Box4A[50],cmap="Greys_r")
plt.subplot(132)
plt.imshow(Box4B[50], cmap="Greys_r")
plt.subplot(133)
plt.imshow(Box4A[50]-Box4B[50], cmap="Greys_r")
plt.show()


