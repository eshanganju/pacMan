#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:01:47 2019

@author: eg
"""

# Imports: 
import numpy as np
from skimage import restoration
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
import spam.datasets as sdata
import spam.kalisphera as skali
import math
import scipy

# %% Making spheres:
pixelSize = 30.e-6                                                      # m/pixel


blurSTD = 0.8                                                           # The standard deviation of the image blurring to be applied first
noiseSTD = 0.03                                                         # The standard deviation of the random noise to be added to each voxel


boxSizeDEM, centres, radii = sdata.loadDEMboxsizeCentreRadius()         # Loading data from SPAM database - centres, radii and box size (length units)
rMax = np.amax(radii)                                                   # Get maximum radius to pad our image (periodic boundaries...)
boxSize = boxSizeDEM + 3 * rMax                                         # Make the box such that there are no particle cutting the edges
centres[:, :] = centres[:, :] + 1.5 * rMax                              # Move the positions to the new center of the image
boxSize = int(math.ceil(np.max(boxSize[:]) / pixelSize))
centres = centres / pixelSize
radii = radii / pixelSize
Box = np.zeros((boxSize, boxSize, boxSize), dtype="<f8")                # Box holding the 3D CT data, Datatype is a 64-bit floating-point number, "<" not sure
skali.makeSphere(Box, centres, radii)                                   # Kalisphera [1] used to genetrate the spheres - Accounting for PVE. particle = 1, void = 0
Box1 = Box
Box1[np.where(Box1 > 1.0)] = 1.0                                        # Checking for outliers voxels x|1<x; DEM uses overlap to compute forces hence vox > 1
Box1[np.where(Box1 < 0.0)] = 0.0                                        # Checking for outliers voxels x|0>x; Why this will be needed is unclear
Box2 = Box1 * 0.5                                                       # All voxles are squeezed from 0-1 to 0-0.5 
Box2 = Box2 + 0.25                                                      # All voxels are shifted from 0-0.5 to 0.25-0.75


Box3 = scipy.ndimage.filters.gaussian_filter(Box2, sigma=blurSTD)       # Applies a gaussian blur. Good example in [3]
Box4 = np.random.normal(Box3, scale=noiseSTD)                           # Applies a random noise to the data [4]


oriClean = Box2[50]
ori50 = Box4[50]


# %% Filters
ori50_cham = restoration.denoise_tv_chambolle(ori50)
ori50_median = median_filter(ori50, size=3)
ori50_nlm3 = restoration.denoise_nl_means(ori50, patch_size=3)
ori50_nlm5 = restoration.denoise_nl_means(ori50, patch_size=5)


# Plotting
ax1=plt.subplot(231)
ax1 = plt.imshow(oriClean, cmap='Greys_r')
plt.title('Original with PVE')


ax2=plt.subplot(232)
ax2.imshow(ori50, cmap='Greys_r')
plt.title('PVE+noise+blur')

ax3=plt.subplot(233)
ax3.imshow(ori50_cham, cmap='Greys_r')
plt.title('Denoised - Chambolle')


ax4=plt.subplot(234)
ax4.imshow(ori50_median, cmap='Greys_r')
plt.title('Denoised - median')


ax5=plt.subplot(235)
ax5.imshow(ori50_nlm3, cmap='Greys_r')
plt.title('Denoised - nlm3')


ax6=plt.subplot(236)
ax6.imshow(ori50_nlm5, cmap='Greys_r')
plt.title('Denoised - nlm5')

plt.savefig('plts.tiff',dpi=600)
plt.show()
plt.close()


# %% Histograms
num = (oriClean.shape[0]*oriClean.shape[0])
lis_oriClean = oriClean.reshape((num,1))
lis_ori50 = ori50.reshape((num,1))
lis_ori50_cham = ori50_cham.reshape((num,1))
lis_ori50_median = ori50_median.reshape((num,1))
lis_ori50_nlm3 = ori50_nlm3.reshape((num,1))
lis_ori50_nlm5 = ori50_nlm5.reshape((num,1))

his_oriClean = np.histogram(lis_oriClean, bins = 1000, range=(0,1))
his_ori50 = np.histogram(lis_ori50, bins = 1000, range=(0,1))
his_ori50_cham = np.histogram(lis_ori50_cham, bins = 1000, range=(0,1))
his_ori50_median = np.histogram(lis_ori50_median, bins = 1000, range=(0,1))
his_ori50_nlm3 = np.histogram(lis_ori50_nlm3, bins = 1000, range=(0,1))
his_ori50_nlm5 = np.histogram(lis_ori50_nlm5, bins = 1000, range=(0,1))
his_x = np.arange(0.0005,1,0.001)

ax7 = plt.subplot(111)
ax7.plot(his_x, his_ori50[0])
ax7.plot(his_x, his_ori50_cham[0])
ax7.plot(his_x, his_ori50_median[0])
ax7.plot(his_x, his_ori50_nlm3[0])
ax7.plot(his_x, his_ori50_nlm5[0])
plt.legend(['Noisy','Cham','Med','nlm3','nlm5'])
plt.xlim((0,1))
plt.ylim((0,600))

plt.savefig('hist.tiff',dpi=600)
plt.show()
plt.close()
