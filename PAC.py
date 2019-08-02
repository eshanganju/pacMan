# -*- coding: utf-8 -*-

'''
PAC - Particle Analysis Code

Analyze static CT scans for particle properties

Features:
    1. Filter tiff stack
    2. Binarize tiff stack
    3. Segment tiff stack
    4. Determine particle size parameters
    5. Determine particle morphology parameters
    6. Determine particle contacts and fabric tensor
    7. Determine particle orientation

Classes:
    1. Particle
    2. Aggregate
---
    3. Filter
    4. Segment
    5. Measure
    6. Visualize

Naming conventions:
    Class names: concatenated words each starting with upper case
    Objects, ivars, methods: concatenated words, first word all lower case, subsequent words starting with upper case (camel case)

Reference [All reference will be cited in the method or class]: 
    []
    []

EG 2019
'''

# %% Importing libraries

# Routine py3 libraries
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.ndimage

# Spam datasets
import spam.datasets as sdata
import spam.kalisphera as skali

# Import class libraries: 
import Particle             # Used to store informaiton about particle
import Aggregate            # Used to store information about aggregate
import Filter               # Contains methods that filter image data
import Segment              # Contains methods that segment particles
import Measure              # Contains methods that particle and aggregate properties
import Visualize            # Contains methods to plot data

# %% DATA option #1 - Kalisphera spheres

# m/pixel
pixelSize = 30.e-6

# The standard deviation of the image blurring to be applied first
blurSTD = 0.8

# The standard deviation of the random noise to be added to each voxel
noiseSTD = 0.03
boxSizeDEM, centres, radii = sdata.loadDEMboxsizeCentreRadius()

# get maximum radius to pad our image (periodic boundaries...)
rMax = np.amax(radii)
boxSize = boxSizeDEM + 3 * rMax

# move the positions to the new center of the image
centres[:, :] = centres[:, :] + 1.5 * rMax

# turn the mm measures into pixels
boxSize = int(math.ceil(np.max(boxSize[:]) / pixelSize))
centres = centres / pixelSize
radii = radii / pixelSize

Box = np.zeros((boxSize, boxSize, boxSize), dtype="<f8")
skali.makeSphere(Box, centres, radii)

# Kalisphera image - blk and white
Box[np.where(Box > 1.0)] = 1.0
Box[np.where(Box < 0.0)] = 0.0
plt.figure()
plt.imshow(Box[Box.shape[0] // 2], vmin=0, vmax=1, cmap='Greys_r')
plt.title("Kalisphera")
plt.show()

# Re-setting image
Box2 = Box * 0.5
Box2 = Box2 + 0.25
plt.figure()
plt.imshow(Box2[Box2.shape[0] // 2], vmin=0, vmax=1, cmap='Greys_r')
plt.title("Kalisphera rescaled [0.25, 0.75]")
plt.show()

# Gaussian filtering
Box3 = scipy.ndimage.filters.gaussian_filter(Box2, sigma=blurSTD)
plt.figure()
plt.imshow(Box3[Box3.shape[0] // 2], vmin=0, vmax=1, cmap='Greys_r')
plt.title("Kalisphera blurred")
plt.show()

# Noise
Box4 = np.random.normal(Box3, scale=noiseSTD)
plt.figure()
plt.imshow(Box4[Box4.shape[0] // 2], vmin=0, vmax=1, cmap='Greys_r')
plt.title("Kalisphera with noise")
plt.show()

# %% DATA option #2 - CT data in aggregate object


# %% Filtering image data in aggregate

# Create filter object
flt = Filter()


# %% Segment CT data

# Creater segmet object
seg = Segment()


# %% Measure data [size, morphology, fabric]

# Create measure object
mes = Measure()

# Size

# Morphology

# Fabric


#%% Visualize [Size, porphology and fabric]

# Creater a visulaize object
viz = Visualize()

# Size

# Morphology

# Fabric



