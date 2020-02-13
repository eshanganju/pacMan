# Imports: 
import numpy as np
import spam.datasets as sdata
import spam.kalisphera as skali
import math
import skimage.external.tifffile as tiffy
import scipy

dem = False
single = False
double = True


if dem == True:
    # DEM Data
    pixelSize = 10e-6       # m/pixel
    blurSTD = 0.8           # std blur
    noiseSTD = 0.03         # Noise

    boxSizeDEM, centres, radii = sdata.loadDEMboxsizeCentreRadius()
    rMax = np.amax(radii)
    boxSize = boxSizeDEM + 3 * rMax
    centres[:, :] = centres[:, :] + 1.5 * rMax

    boxSize = int(math.ceil(np.max(boxSize[:]) / pixelSize))
    centres = centres / pixelSize
    radii = radii / pixelSize
    Box = np.zeros((boxSize, boxSize, boxSize), dtype=np.uint32)

    skali.makeSphere(Box, centres, radii)
    Box1 = Box
    Box1[np.where(Box1 > 1.0)] = 1.0
    Box1[np.where(Box1 < 0.0)] = 0.0
    Box2 = Box1 * 0.5
    Box2 = Box2 + 0.25
    Box3 = scipy.ndimage.filters.gaussian_filter(Box2, sigma=blurSTD)
    Box4 = np.random.normal(Box3, scale=noiseSTD) 
    tiffy.imsave('Box4B.tiff',Box4.astype('uint32'))

if single == True:
    # Single particle
    boxSize = [100,100,100]
    centres = [49,49,49]
    radii=[25]
    Box = np.zeros((boxSize, boxSize, boxSize), dtype=np.uint32)
    skali.makeSphere(Box, centres, radii)
    Box1 = Box
    Box1[np.where(Box1 > 1.0)] = 1.0
    Box1[np.where(Box1 < 0.0)] = 0.0
    Box2 = Box1 * 0.5
    Box2 = Box2 + 0.25
    Box3 = scipy.ndimage.filters.gaussian_filter(Box2, sigma=blurSTD)
    Box4 = np.random.normal(Box3, scale=noiseSTD) 
    tiffy.imsave('single.tiff',Box4.astype('uint32'))


if single == True:
    # Single particle
    boxSize = [200,100,100]
    centres = [[74,49,49],[124,49,49]]
    radii=[25,25]
    Box = np.zeros((boxSize, boxSize, boxSize), dtype=np.uint32)
    skali.makeSphere(Box, centres, radii)
    Box1 = Box
    Box1[np.where(Box1 > 1.0)] = 1.0
    Box1[np.where(Box1 < 0.0)] = 0.0
    Box2 = Box1 * 0.5
    Box2 = Box2 + 0.25
    Box3 = scipy.ndimage.filters.gaussian_filter(Box2, sigma=blurSTD)
    Box4 = np.random.normal(Box3, scale=noiseSTD) 
    tiffy.imsave('double.tiff',Box4.astype('uint32'))

