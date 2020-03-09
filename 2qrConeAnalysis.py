'''
I want o extract a subvolume knowing the size and location of the top left corner of the subregion
I want to get binary and segmented images

'''

import numpy as np
import matplotlib.pyplot as plt

from classes import Aggregate
from classes import Reader
from classes import Filter
from classes import Segment
import spam.label as slab
import skimage.external.tifffile as tf
#from classes import Measure

outputDirectory = '/home/eg/codes/pacOutput/'
tiffFileLocation = '/home/eg/codes/pacInput/2QR_90_tip/'

sampleName = '2qr-90kPa-tip'
centerSliceOfdata = 513
subregionTopLeftY = 711
subregionTopLeftX = 171

calib = 0.0122996 # mm/px
d50 = 0.72 # mm
sizeOfSubregionCube = 5*d50
gliMax = 2**16-1

centerZ = centerSliceOfdata
centerY = int(subregionTopLeftY + sizeOfSubregionCube/calib//2)
centerX = int(subregionTopLeftX + sizeOfSubregionCube/calib//2)

unfilteredData = Reader.readTiffFileSequence(tiffFileLocation,centerZ,centerY,centerX,sizeOfSubregionCube,calib)
#plt.imshow(unfilteredData[unfilteredData.shape[0]//2],cmap='gray')
#plt.show()

#filteredData = Filter.filterDenoiseNlm(sampleName,unfilteredData,gliMax,outputDirectory)
#plt.imshow(filteredData[filteredData.shape[0]//2])
#plt.show()

binaryData, otsuThreshold = Segment.binarizeAccordingToOtsu(unfilteredData,outputDirectory,sampleName)
#plt.imshow(binaryData[binaryData.shape[0]//2],cmap='gray')
#plt.show()

labeledData = slab.watershed(binaryData)

labelledMapName = outputDirectory + sampleName + '-spamLabelledMap.tiff'
binaryMapName = outputDirectory + sampleName + '-binMap.tiff'
unfilteredMapName = outputDirectory + sampleName + '-gliMap.tiff'

tf.imsave(labelledMapName, labeledData.astype('uint32'))
tf.imsave(binaryMapName, binaryData.astype('uint32'))
tf.imsave(unfilteredMapName, unfilteredData.astype('uint32'))

