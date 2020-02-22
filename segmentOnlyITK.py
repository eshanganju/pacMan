import spam.label as slab
import skimage.external.tifffile as tf
from classes import Segment

s = Segment.Segment()

knownVoidRatio = 0.4983
outputLocation = '/home/eg/codes/pacOutput/'
sampleName = 'otc-1500n'

fgli = tf.imread('/home/eg/codes/pacOutput/otc-1500n-fgliMap.tiff')
binMap = s.binarizeAccordingToDensity(fgli,knownVoidRatio,outputLocation,sampleName)

labMap = slab.watershed(binMap) 

tf.imsave('/home/eg/codes/pacOutput/specialSPAM/otc-1500n-labMapITKWatershed.tiff',labMap)
