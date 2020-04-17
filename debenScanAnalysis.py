'''
The objective of this code is to assemble the codes in pac directory to:
    1. Read the data files for the deben scans and extract volume
    2. Segment the data (binarize and WS)
    3. Analyze CT data for
        3.1 Correct REV size
        3.2 Particle size for REV size
        3.3 Particle morphology
        3.4 Relative breakage
        3.5 Fabric tensor
    4. Distribution of rel. breakage and fabric tensor
'''
# Importing files
from pac import Reader                      # Reads files, cuts into smaller pieces
from pac import Filter                      # Filters files using NLM filter
from pac import Segment                     # Binarizatio and WS
from pac import Measure                     # Calculates particle size, mprphology, contact, breakage
import matplotlib.pyplot as plt             # matplotlib
import skimage.external.tifffile as tf      # scikit-image
import numpy as np                          # numpy

'''
Binarization according to density measurement
Segmentation according to ITK and auto correction
Particle size values and gradation
Relative breakage according to Einav
Contact according to ITK and RW
'''

# Data location:
inputFolderLocation = '/home/eg/codes/pacInput/OTC-0N/'
outputFolderLocation = '/home/eg/codes/pacOutput/OTC-0N/'

# Data details:
dataName = 'otc-0N'
measuredVoidRatio = 0.541   # Void ratio measured from 1D compression experiment
d50 = 0.72                  # D50 in mm - original gradation
calib = 0.01192             # calibration from CT mm/voxel
zCenter = 513               # Voxel units - center of slice
yCenter = 457               # Voxel units - vertical center
xCenter = 507               # Voxel units - horizontal center
edgeLengthMin = 5           # Min edge length in D50s
edgeLengthMax = 6           # Max edge length in D50s


# Reading and cropping the data file
for edgeLength in range( edgeLengthMin , edgeLengthMax ):

    gliMap = Reader.readTiffFileSequence( inputFolderLocation, zCenter, yCenter, xCenter, edgeLength, calib)
    labMap = Segment.obtainLabelledMapUsingITKWS( gliMap , measuredVoidRatio )
    corLabMap = Segment.fixErrorsInSegmentation( labMap )
    neCorLabMap = Segment.removeEdgeLabels( corLabMap )
    gsd1, gsd2, gsd3, gsd4 = Measure.gsd( neCorLabMap )

    '''
    relativeBreakage
    Contactnormals and fabric
    '''

    # Save files as csv
    np.savetxt((outputFolderLocation+ str(edgeLength) +'D50-gsd1.csv'), gsd1, delimiter=',') # Eqsp
    np.savetxt((outputFolderLocation+ str(edgeLength) +'D50-gsd2.csv'), gsd2, delimiter=',') # CA max
    np.savetxt((outputFolderLocation+ str(edgeLength) +'D50-gsd3.csv'), gsd3, delimiter=',') # CA med
    np.savetxt((outputFolderLocation+ str(edgeLength) +'D50-gsd4.csv'), gsd4, delimiter=',') # CA min

    # Save the 3D maps as tiff
    gliName = outputFolderLocation + 'gli.tiff'
    labName = outputFolderLocation + 'lab.tiff'
    corLabName = outputFolderLocation + 'labCor.tiff'
    neCorLabName = outputFolderLocation + 'neCorLab.tiff'

    tf.imsave(gliName,gliMap.astype('uint32'))
    tf.imsave(labName,labMap.astype('uint32'))
    tf.imsave(corLabName,corLabMap.astype('uint32'))
    tf.imsave(neCorLabName,neCorLabMap.astype('uint32'))

