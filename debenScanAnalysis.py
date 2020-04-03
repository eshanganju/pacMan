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

#print('Maithilee chutiya hai')
#print('Haan Haan, barobar')

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
    gsd1, gsd2, gsd3, gsd4 = Measure.gsd(neCorLabMap)

    '''
    DO:
        2. Get relative breakage for different size param
        3. Get contact vectors, fabric tensor and plots for contact
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


'''
# Important to make the starting and ending points the same and x and log x before any changes
x1 = np.arange(0.1,1001,100)
logX1 = np.log10(x1)
y1 = x1**2

x2 = np.arange(0.1,1001,100)
logX2 = np.log10(x2)
y2 = np.array([0,0,0,0,0,0,0,y1.max()//4,y1.max()//2,y1.max()*3//4,y1.max()])

# Area between curves:
bins = 1000
incrLogX = ( logX1.max() - logX1.min() ) / bins
logXInterp = np.arange(logX1.min(),logX1.max(),incrLogX)
y1Interp = np.zeros_like(logXInterp)
y2Interp = np.zeros_like(logXInterp)
area = np.zeros_like(logXInterp)

for i in range(0,bins):
    y1Interp[i] = np.interp(logXInterp[i], logX1, y1 )
    y2Interp[i] = np.interp(logXInterp[i], logX2, y2 )

    if i == 0:
        area[i] = (y1Interp[i] - y2Interp[i])*incrLogX/2
    elif i == bins - 1:
        area[i] = (y1Interp[i] - y2Interp[i])*incrLogX/2
    else:
        area[i] = (y1Interp[i] - y2Interp[i])*incrLogX

totalArea = np.sum(area)
'''


