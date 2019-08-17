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

f = open("logFileAlphaFilter.txt","w+")

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

oriClean = Box2[50]
ori50 = Box4[50]

ori50_cham = restoration.denoise_tv_chambolle(ori50)
ori50_median = median_filter(ori50, size=3)


sigma_est = np.mean(restoration.estimate_sigma(ori50))
patch = dict(patch_size=3,
             patch_distance=25)
ori50_nlm3_patch = restoration.denoise_nl_means(ori50,
                                                h=0.8*sigma_est,
                                                sigma=sigma_est,
                                                fast_mode=False,
                                                **patch)
ori50_nlm3 = restoration.denoise_nl_means(ori50, patch_size=3)
ori50_nlm5 = restoration.denoise_nl_means(ori50, patch_size=5)

ax1=plt.subplot(231)
ax1 = plt.imshow(oriClean, cmap='Greys_r')
plt.title('Original with PVE')

ax2=plt.subplot(232)
ax2.imshow(ori50, cmap='Greys_r')
plt.title('PVE+noise+blur')

ax3=plt.subplot(233)
ax3.imshow(ori50_nlm3_patch, cmap='Greys_r')
plt.title('Denoised - NLM3-Patch')


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
plt.xticks(np.arange(0,1.25,step=0.25))
plt.savefig('hist.tiff',dpi=600)

Clean_nlm3_diff = (oriClean - ori50_nlm3)
plt.figure()
plt.imshow(Clean_nlm3_diff, cmap = 'PiYG')
plt.colorbar()
plt.title('different in pixel values (positive - underdetection)')
plt.savefig('clean-nlm3diff.tiff', dpi=600)

meanOfDifference = np.mean(Clean_nlm3_diff)


f.write("meanDifference = %f\n" % meanOfDifference)
f.write("sigma image = %f\n" % sigma_est)
f.close()

print("###########################")
print("Files made:")
print("\t1. hist.tiff")
print("\t2. plts.tiff")
print("\t3. clean-nlm3diff.tiff")
print("\t4. logFileAlphaFilter.txt")
print("###########################")
