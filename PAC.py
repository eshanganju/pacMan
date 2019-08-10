# -*- coding: utf-8 -*-

'''
Outline: 
PAC - Particle Analysis Code I wanted to call it PACman, but cant find any way of using those 3 letters
Contact me ganju.eshan@gmail.com if you do. This code analyzes static CT scans for particle properties 
(does not track particle)

---
Features:
    1. Filter tiff stack of CT data
    2. Binarize tiff stack of CT data
    3. Segment tiff stack of CT data
    4. Determine particle size parameters of segmented data
    5. Determine particle morphology parameters of segmented data
    6. Determine particle contacts and fabric tensor of segmented data
    7. Determine particle and void orientation

---
Classes:
    1. Particle
    2. Aggregate
---
    3. Filter
    4. Segment
    5. Measure
    6. Visualize

---
Naming conventions:
    Class names: concatenated words each starting with upper case
    Objects, ivars, methods: concatenated words, first word all lower case, subsequent words starting with upper case

---
References [All reference will be cited in the method or class]:
    [1] Tengattini, Alessandro, and Edward Andò. 2015. "Kalisphera: An Analytical Tool to Reproduce the Partial Volume Effect of Spheres Imaged in 3D."
    Measurement Science and Technology 26 (9).https://doi.org/10.1088/0957-0233/26/9/095606.
    [2] Kalisphera example: https://ttk.gricad-pages.univ-grenoble-alpes.fr/spam/spam_examples/kalisphera/plot_generateAssembly.html#sphx-glr-spam-examples-kalisphera-plot-generateassembly-py
    [3] Gaussian blur example: https://stackoverflow.com/questions/25216382/gaussian-filter-in-scipy
    [4] Random noise: https://docs.scipy.org/doc/numpy-1.14.2/reference/generated/numpy.random.normal.html

EG 2019
'''

# %% Importing libraries

# Python libraries
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.ndimage

# Spam libraries
import spam.datasets as sdata
import spam.kalisphera as skali

# Contructed libraries
import Aggregate            # Used to store information about aggregate
import Filter               # Contains methods that filter image data
import Segment              # Contains methods that segment particles
import Measure              # Contains methods that particle and aggregate properties
import LemmeC               # Contains methods to plot data
import Writer               # Writer - makes files that can later be used

#---#---#---#
# TODO: Minimize unnesary commenting

#%% Activate tools

cleanUp = Filter.Filter()       # Filter tool
segment = Segment.Segment()     # Segment tool
measure = Measure.Measure()     # Measure tool
lemmeC = LemmeC.LemmeC()        # Visualization tool
writer = Writer.Writer()        # Writer to store analysis data

#---#---#---#
# TODO: Minimize unnesary commenting

# %% DATA option #1 - Kalisphera spheres

'''
Use this segment if testing functionality of code
Use next segment if using collected CT data
---
This section uses data from spam.datasets (sdata).
Data is in the form of DEM spheres: The centeres, radius and box size are taken from loadDEMboxsizeCentreRadius
The data is used to create touching spheres using kalisphera tool developed by [1], which acounts for PVE
The generated data is passed through a gaussian and a noise filter to make it similar to data colllected from X-Ray CT
The example is slightly modified after [2]
'''

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

aggregate = Aggregate.Aggregate(Box4, pixelSize)                        # Box 4 dataset is passed (PVE'd DEM data; blurred and filtered)

measure.measureBenchMarkSizeAndNormal(aggregate,radii,centres)          # Generates grain size distribution and contact normals


# %% DATA option #2 - CT data in aggregate object
'''
In this section we upload data from CT scans 
This will need some pre-processing to filter the images
keep greylevel data between 0 and 1
Visualization wil be very useful in this case. 
'''

#---#---#---#
# TODO: Add visualization at this stage

# %% Segment CT data

segment.greyLevelHistogram(aggregate)                   # Generate the grey level histogram and store in aggregate as "greyLevelHistogram"
segment.binarizeOtsu(aggregate)                         # Binarization using OTSU and store threshold in aggregate as "globalThreshold"
segment.euclidDistanceMap(aggregate)                    # Create EDM and store in aggregate as "euclidDistanceMap"
segment.localMaximaMarkers(aggregate)                   # Locate peaks of EDM, store location of each marker (label) as "markers"
segment.topoWatershed(aggregate)                        # Topological watershed and store in aggregate as "labelledMap" create a "particleList"

#---#---#---#
# TODO: Number of particles keeps changing - check segment.localMaximaMarkers
# TODO: Correction of oversegmentation - long contact edge detection
# TODO: Implement random Walker segmentation, manual checks

# %% Measure data [size, morphology, fabric]

measure.measureParticleSizeDistribution(aggregate)      # Size
measure.measureMorphology(aggregate)                    # Morphology
measure.measureContactNormalsSpam(aggregate)            # Fabric from SPAM libraries

#---#---#---#

#%% Visualize [Size, morphology and fabric]

# Creater a visulaize object

# Size

# Morphology

# Fabric



