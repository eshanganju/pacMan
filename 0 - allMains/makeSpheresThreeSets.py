# Imports: 
import numpy as np
import spam.datasets as sdata
import spam.kalisphera as skali
import math
import skimage.external.tifffile as tiffy
import scipy
import time
import matplotlib.pyplot as plt

dem = False
single = False
double = True


if dem == True:
    start = time.time()

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
    Box = np.zeros((boxSize, boxSize, boxSize), dtype='<f8')

    skali.makeSphere(Box, centres, radii)
    Box1b = Box
    Box1b[np.where(Box1b > 1.0)] = 1.0
    Box1b[np.where(Box1b < 0.0)] = 0.0
    Box2b = Box1b * 0.5
    Box2b = Box2b + 0.25
    Box3b = np.round(Box2b*((2**32)-1))
    Box4b = scipy.ndimage.filters.gaussian_filter(Box3b, sigma=blurSTD)
    Box5b = np.random.normal(Box4b, scale=noiseSTD) 
    Box7b = Box5b.astype('uint32')
    tiffy.imsave('NEW-DEM.tiff',Box5b.astype(np.uint32))
    end = time.time()
    timeTaken = np.round((end-start),2)
    print("Time taken for DEM = %f" % timeTaken)

if single == True:
    start = time.time()
    blurSTD = 0.8
    # Single particle
    boxSize = 100
    centres = [49,49,49]
    radii=[25]
    Box = np.zeros((boxSize, boxSize, boxSize), dtype='<f8')

    skali.makeSphere(Box, centres, radii)
    Box1b = Box
    Box1b[np.where(Box1b > 1.0)] = 1.0
    Box1b[np.where(Box1b < 0.0)] = 0.0
    Box2b = Box1b * 0.5
    Box2b = Box2b + 0.25
    Box3b = np.round(Box2b*((2**32)-1))
    Box4b = scipy.ndimage.filters.gaussian_filter(Box3b, sigma=blurSTD)
    Box5b = np.random.normal(Box4b, scale=noiseSTD) 
    Box7b = Box5b.astype('uint32')
    tiffy.imsave('NEW-Single.tiff',Box5b.astype(np.uint32))
    end = time.time()
    timeTaken = np.round((end-start),2)
    print("Time taken for Single = %f" % timeTaken)

if double == True:
    start = time.time()
    blurSTD = 0.8

    # Single particle
    boxSizeXX = 100
    boxSizeYY = 100
    boxSizeZZ = 200
    centres = [[74,49,49],[124,49,49]]
    radii=[27,27]
    Box = np.zeros((boxSizeZZ, boxSizeYY, boxSizeXX), dtype='<f8')

    skali.makeSphere(Box, centres, radii)
    Box1b = Box
    Box1b[np.where(Box1b > 1.0)] = 1.0
    Box1b[np.where(Box1b < 0.0)] = 0.0
    Box2b = Box1b * 0.5
    Box2b = Box2b + 0.25
    Box3b = np.round(Box2b*((2**32)-1))
    Box4b = scipy.ndimage.filters.gaussian_filter(Box3b, sigma=blurSTD)
    Box5b = np.random.normal(Box4b, scale=blurSTD) 
    Box7b = Box5b.astype('uint32')
    tiffy.imsave('NEW-Double.tiff',Box5b.astype(np.uint32))
    end = time.time()
    timeTaken = np.round((end-start),2)
    print("Time taken for Double = %f" % timeTaken)

