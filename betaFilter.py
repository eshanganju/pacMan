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

sigma_est = np.mean(restoration.estimate_sigma(ori50, multichannel=True))

nlm_3ps_11pd_02h = restoration.denoise_nl_means(ori50, 
                                                patch_size=3, 
                                                patch_distance=11, 
                                                h=0.02, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)

nlm_5ps_11pd_02h = restoration.denoise_nl_means(ori50, 
                                                patch_size=5, 
                                                patch_distance=11, 
                                                h=0.02, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)

nlm_7ps_11pd_02h = restoration.denoise_nl_means(ori50, 
                                                patch_size=7, 
                                                patch_distance=11, 
                                                h=0.02, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)

ax1=plt.subplot(131)
ax1.imshow(nlm_3ps_11pd_02h, cmap='Greys_r')
plt.title('NLM3ps11pd0.02h')

ax2=plt.subplot(132)
ax2.imshow(nlm_5ps_11pd_02h, cmap='Greys_r')
plt.title('NLM5ps11pd0.02h')

ax3=plt.subplot(133)
ax3.imshow(nlm_7ps_11pd_02h, cmap='Greys_r')
plt.title('NLM7ps11pd0.02h')

plt.savefig('plts_nlm.tiff',dpi=600) 
plt.show()

nlm_5ps_09pd_02h = restoration.denoise_nl_means(ori50, 
                                                patch_size=5, 
                                                patch_distance=9, 
                                                h=0.02, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)

nlm_5ps_11pd_02h = restoration.denoise_nl_means(ori50, 
                                                patch_size=5, 
                                                patch_distance=11, 
                                                h=0.02, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)

nlm_5ps_13pd_02h = restoration.denoise_nl_means(ori50, 
                                                patch_size=7, 
                                                patch_distance=13, 
                                                h=0.02, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)

ax1=plt.subplot(131)
ax1.imshow(nlm_5ps_09pd_02h, cmap='Greys_r')
plt.title('NLM5ps09pd0.02h')

ax2=plt.subplot(132)
ax2.imshow(nlm_5ps_11pd_02h, cmap='Greys_r')
plt.title('NLM5ps11pd0.02h')

ax3=plt.subplot(133)
ax3.imshow(nlm_5ps_13pd_02h, cmap='Greys_r')
plt.title('NLM5ps13pd0.02h')

plt.savefig('plts_nlm_2.tiff',dpi=600) 
plt.show()

