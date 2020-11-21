import spam.label as slab
import skimage.external.tifffile as tf
from classes import Segment


knownVoidRatio = 0.5405
outputLocation = '/home/eg/codes/pacOutput/'
sampleName = 'otc-0n-small'

fgli = tf.imread('/home/eg/codes/pacOutput/otc-0n-small-fgliMap.tiff')
binMap = Segment.binarizeAccordingToDensity(fgli,knownVoidRatio,outputLocation,sampleName)

labMap = slab.watershed(binMap) 

tf.imsave('/home/eg/codes/pacOutput/specialSPAM/otc-0n-small-labMapITKWatershed.tiff',labMap)
