# Importing files
from pac import Reader                      # Reads files, cuts into smaller pieces
from pac import Filter                      # Filters files using NLM filter
from pac import Segment                     # Binarizatio and WS
from pac import Measure                     # Calculates particle size, mprphology, contact, breakage
from pac import Plot                        # Plotting functions
import time                                 # The fourth dimension
import matplotlib.pyplot as plt             # matplotlib
import skimage.external.tifffile as tf      # scikit-image
import numpy as np    


inputFolderLocation = '/home/eg/codes/pacInput/OGF-1500N/'
outputFolderLocation = '/home/eg/codes/pacOutput/OGF-1500N/'
originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogf30MPa.csv' # Original GSD location


# Data details 0N:
dataName = 'ogf-1500N'
measuredVoidRatioSample = 0.565                                     # Void ratio measured from 1D compression experiment
d50 = 0.62                                                          # D50 in mm - original gradation
cal = 0.01193                                                       # calibration from CT mm/voxel
zCenter = 513                                                       # Voxel units - center of slice
yCenter = 458                                                       # Voxel units - vertical center
xCenter = 492                                                       # Voxel units - horizontal center
originalGSD = np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

edgeLength = 5 

gliMap = Reader.readTiffFileSequence( inputFolderLocation, zCenter, yCenter, xCenter, edgeLength, cal, invertImageData=False)

binMap, edMap, edPeaksMap, labMap = Segment.obtainLabelledMapUsingITKWS( gliMap , measuredVoidRatio=measuredVoidRatioSample , outputLocation=outputFolderLocation )
correctedLabMap = Segment.fixErrorsInSegmentation( labMap , pad=2)

tf.imsave('gliMap.tiff',gliMap.astype('uint32'))
tf.imsave('binMap.tiff',binMap.astype('uint32'))
tf.imsave('edMap.tiff',edMap.astype('uint32'))
tf.imsave('edPeaksMap.tiff',edPeaksMap.astype('uint32'))
tf.imsave('labMap.tiff',labMap.astype('uint32'))
tf.imsave('corLabMap.tiff',correctedLabMap.astype('uint32'))
