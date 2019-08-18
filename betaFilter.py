# Imports: 
import numpy as np
from skimage import restoration
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
import spam.datasets as sdata
import spam.kalisphera as skali
import math
import scipy

f = open("logFileBetaFilter.txt","w+")

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

axx1 = plt.subplot(121)
axx1.imshow(oriClean, cmap="Greys_r")
plt.title("Clean + PVE")
axx2 = plt.subplot(122)
axx2.imshow(ori50, cmap="Greys_r")
plt.title("PVE+Blur+Noise")
plt.savefig('pltsEffectOfNoise.tiff',dpi=600)


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

nlm_5ps_09pd_02h = restoration.denoise_nl_means(ori50, 
                                                patch_size=5, 
                                                patch_distance=9, 
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

nlm_5ps_11pd_05h = restoration.denoise_nl_means(ori50, 
                                                patch_size=5, 
                                                patch_distance=11, 
                                                h=0.05, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)

nlm_5ps_11pd_10h = restoration.denoise_nl_means(ori50, 
                                                patch_size=5, 
                                                patch_distance=11,
                                                h=0.1, 
                                                multichannel=False, 
                                                fast_mode=True, 
                                                sigma=sigma_est)
noiseRemoved = ori50 - nlm_5ps_11pd_05h

plt.figure()
ax0=plt.imshow(noiseRemoved, cmap='Greys_r')
plt.title("Noise removed")
plt.colorbar(ax0)
plt.savefig('pltsNoise.tiff',dpi=600)

ax1=plt.subplot(331)
ax1.imshow(nlm_3ps_11pd_02h, cmap='Greys_r')
plt.title('NLM3ps11pd0.02h')

ax2=plt.subplot(332)
ax2.imshow(nlm_5ps_11pd_02h, cmap='Greys_r')
plt.title('NLM5ps11pd0.02h')

ax3=plt.subplot(333)
ax3.imshow(nlm_7ps_11pd_02h, cmap='Greys_r')
plt.title('NLM7ps11pd0.02h')

ax4=plt.subplot(334)
ax4.imshow(nlm_5ps_09pd_02h, cmap='Greys_r')
plt.title('NLM5ps09pd0.02h')

ax5=plt.subplot(335)
ax5.imshow(nlm_5ps_11pd_02h, cmap='Greys_r')
plt.title('NLM5ps11pd0.02h')

ax6=plt.subplot(336)
ax6.imshow(nlm_5ps_13pd_02h, cmap='Greys_r')
plt.title('NLM5ps13pd0.02h')

ax7=plt.subplot(337)
ax7.imshow(nlm_5ps_11pd_02h, cmap='Greys_r')
plt.title('NLM5ps11pd0.02h')

ax8=plt.subplot(338)
ax8.imshow(nlm_5ps_11pd_05h, cmap='Greys_r')
plt.title('NLM5ps11pd0.05h')

ax9=plt.subplot(339)
ax9.imshow(nlm_5ps_11pd_10h, cmap='Greys_r')
plt.title('NLM5ps11pd0.10h')

plt.savefig('pltsNLM.tiff',dpi=600) 

# histogram
num = (oriClean.shape[0]*oriClean.shape[0])
listOri50 = ori50.reshape((num,1))
listNlm05ps11pd02h = nlm_5ps_11pd_02h.reshape((num,1))
listNlm05ps11pd05h = nlm_5ps_11pd_05h.reshape((num,1))
listNlm05ps11pd10h = nlm_5ps_11pd_10h.reshape((num,1))
his_x = np.arange(0.0005,1,0.001)

histOri50 = np.histogram(listOri50, bins=1000, range=(0,1))
histNlm05ps11pd02h = np.histogram(listNlm05ps11pd02h, bins=1000, range=(0,1))
histNlm05ps11pd05h = np.histogram(listNlm05ps11pd05h, bins=1000, range=(0,1))
histNlm05ps11pd10h = np.histogram(listNlm05ps11pd10h, bins=1000, range=(0,1))

ax10 = plt.subplot(111)
ax10.plot(his_x, histOri50[0])
ax10.plot(his_x, histNlm05ps11pd02h[0])
ax10.plot(his_x, histNlm05ps11pd05h[0])
ax10.plot(his_x, histNlm05ps11pd10h[0])
plt.legend(['Noisy','histNlm05ps11pd02h','histNlm05ps11pd05h','histNlm05ps11pd10h'])
plt.xlim((0,1))
plt.ylim((0,600))
plt.xticks(np.arange(0,1.25,step=0.25))
plt.savefig('histNLM.tiff',dpi=600)

# Output files:
print("###########################")
print("Files made:")
print("\t1. histNLM.tiff")
print("\t2. pltsNLM.tiff")
print("\t3. pltsNoise.tiff")
print("###########################")
