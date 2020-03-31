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
import matplotlib.pyplot as plt             #
import skimage.external.tifffile as tf      #

# Data location:
inputFolderLocation = '/home/eg/codes/pacInput/OTC-0N/'
outputFolderLocation = '/home/eg/codes/pacOutput/OTC-0N/'

# Data details:
dataName = 'otc-0N'
d50 = 0.72              # D50 in mm
calib = 0.01192         # calibration from CT mm/voxel
zCenter = 513           # Voxel units
yCenter = 457           # Voxel units
xCenter = 507           # Voxel units
edgeLengthMin = 5       # Min edge length in D50s
edgeLengthMax = 7       # Max edge length in D50s


# Reading and cropping the data file
for edgeLength in range( edgeLengthMin , edgeLengthMin + 1 ):

    gliMap = Reader.readTiffFileSequence( inputFolderLocation, zCenter, yCenter, xCenter, edgeLength, calib)
    labMap = Segment.obtainLabelledMapUsingITKWS(gliMap)
    corLabMap = Segment.fixErrorsInSegmentation(labMap, pad = 2)
    neCorLabMap = Segment.removeEdgeLabels(corLabMap)

    '''
     
    '''

    #print('\nLabelled map: ')
    #Segment.countEdgeLabels(labMap)
    #print('\nCorr Lab map: ')
    #Segment.countEdgeLabels(labCorrMap)
    #print('\nNo Edge Corr Lab map: ')
    #Segment.countEdgeLabels(noEdgeLabCorrMap)

    gliName = outputFolderLocation + 'gli.tiff'
    labName = outputFolderLocation + 'lab.tiff'
    corrLabName = outputFolderLocation + 'labCorr.tiff'
    neCorLabName = outputFolderLocation + 'nelabCorr.tiff'

    tf.imsave(gliName,gliMap.astype('uint32'))
    tf.imsave(labName,labMap.astype('uint32'))
    tf.imsave(corrLabName,corLabMap.astype('uint32'))

