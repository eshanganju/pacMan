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

tiffFileLocation = '/home/eg/codes/pacInput/OGF_25_top/'

sampleName = 'ogf-25kPa-top'
centerSliceOfdata = 513
subregionTopLeftY = 345
subregionTopLeftX = 258

calib = 0.0122996 # mm/px
d50 = 0.62 # mm
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

tf.imsave(labelledMapName, labeledData.astype('uint16'))
tf.imsave(binaryMapName, binaryData.astype('uint16'))
tf.imsave(unfilteredMapName, unfilteredData.astype('uint16'))

